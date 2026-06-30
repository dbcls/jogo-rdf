#!/usr/bin/env ruby
# frozen_string_literal: true
#
# jogo2ttl.rb — JoGo haplotype JSON → RDF/Turtle 変換 (rdf-config の model.yaml に準拠)
#
# TtlBuilder は bin/regions2ttl.rb (本番の一括変換) から再利用される。
# 下記 CLI は単一遺伝子の動作確認用 (emit_common: true で共通ノードもインライン出力)。
#
# 使い方 (プロジェクト直下から):
#   ruby bin/jogo2ttl.rb --genename ALDH2
#   ruby bin/jogo2ttl.rb --regionname ALDH2_chr12_111761933_111822532
#   出力先: 省略時 output/<gene_key>.ttl  (--out PATH で指定, "-" で標準出力)
#
# 取得API:
#   https://jogo.csml.org/gene?genename=ALDH2&format=json
#   https://jogo.csml.org/genicregion?regionname=ALDH2_chr12_111761933_111822532&format=json
#
# 将来: 遺伝子リスト取得APIができたら、その一覧を Jogo.fetch_by_genename に回して
#       全遺伝子をループ変換する (末尾の main 参照)。
#
# 注意: 集団参照データ (population.ttl) と オントロジー (ontology.ttl) は静的なので
#       本スクリプトは出力しない。生成 ttl と併せてロードする想定。

require 'json'
require 'net/http'
require 'uri'
require 'optparse'

BASE = 'https://jogo.csml.org'
SUPERPOPS = %w[AFR AMR EAS EUR SAS].freeze

module Jogo
  module_function

  def fetch_json(url)
    res = Net::HTTP.get_response(URI(url))
    raise "HTTP #{res.code} for #{url}" unless res.is_a?(Net::HTTPSuccess)
    JSON.parse(res.body)
  end

  def fetch_by_genename(name)
    fetch_json("#{BASE}/gene?genename=#{URI.encode_www_form_component(name)}&format=json")
  end

  def fetch_by_regionname(name)
    fetch_json("#{BASE}/genicregion?regionname=#{URI.encode_www_form_component(name)}&format=json")
  end

  # Turtle 文字列リテラルのエスケープ
  def lit(s)
    e = s.to_s.gsub(/[\\"\n\r\t]/) do |c|
      { "\\" => '\\\\', '"' => '\\"', "\n" => '\\n', "\r" => '\\r', "\t" => '\\t' }[c]
    end
    %("#{e}")
  end
end

class TtlBuilder
  PREFIXES = {
    'jogo'      => 'http://jogo.csml.org/ontology#',
    'jogo_gene' => 'http://jogo.csml.org/rdf/gene/',
    'jogo_hap'  => 'http://jogo.csml.org/rdf/haplotype/',
    'jogo_var'  => 'http://jogo.csml.org/rdf/variant/',
    'jogo_pop'  => 'http://jogo.csml.org/rdf/population/',
    'hco'       => 'http://identifiers.org/hco/',
    'hgnc'      => 'http://identifiers.org/hgnc/',
    'ensembl'   => 'http://identifiers.org/ensembl/',
    'faldo'     => 'http://biohackathon.org/resource/faldo#',
    'gvo'       => 'http://genome-variation.org/resource#',
    'foaf'      => 'http://xmlns.com/foaf/0.1/',
    'dcterms'   => 'http://purl.org/dc/terms/',
    'skos'      => 'http://www.w3.org/2004/02/skos/core#',
    'rdf'       => 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
    'rdfs'      => 'http://www.w3.org/2000/01/rdf-schema#'
  }.freeze

  # section, class, idカラム, lenカラム, 複合→個別リンク [[predicate, 個別idカラム], ...]
  SUMMARY = [
    ['ahaplotypesummary',    'AHaplotype',    'ahapid',    'alen', []],
    ['chaplotypesummary',    'CHaplotype',    'chapid',    'clen', []],
    ['thaplotypesummary',    'THaplotype',    'thapid',    'tlen', []],
    ['ghaplotypesummary',    'GHaplotype',    'ghapid',    nil,    []],
    ['achaplotypesummary',   'ACHaplotype',   'achapid',   nil,
     [['aHaplotype', 'ahapid'], ['cHaplotype', 'chapid']]],
    ['acthaplotypesummary',  'ACTHaplotype',  'acthapid',  nil,
     [['aHaplotype', 'ahapid'], ['cHaplotype', 'chapid'], ['tHaplotype', 'thapid']]],
    ['actghaplotypesummary', 'ACTGHaplotype', 'actghapid', nil,
     [['aHaplotype', 'ahapid'], ['cHaplotype', 'chapid'], ['tHaplotype', 'thapid'], ['gHaplotype', 'ghapid']]]
  ].freeze

  # 遺伝子 → 各レベル一覧の述語
  GENE_HAP_PRED = {
    'ahaplotypesummary' => 'aHaplotype', 'chaplotypesummary' => 'cHaplotype',
    'thaplotypesummary' => 'tHaplotype', 'ghaplotypesummary' => 'gHaplotype',
    'achaplotypesummary' => 'acHaplotype', 'acthaplotypesummary' => 'actHaplotype',
    'actghaplotypesummary' => 'actgHaplotype'
  }.freeze

  # variant表: section, レベル, hapid列 (c/t/g は複合id, a は個別id)
  VARIANT_TABLES = [
    ['avariants', 'a', 'ahapids'],
    ['cvariants', 'c', 'chapids'],
    ['tvariants', 't', 'thapids'],
    ['gvariants', 'g', 'ghapids']
  ].freeze

  def self.prefix_header
    PREFIXES.map { |p, u| "@prefix #{p}: <#{u}> ." }.join("\n") + "\n\n"
  end

  def self.block(subject, triples)
    "#{subject}\n    #{triples.join(" ;\n    ")} .\n\n"
  end

  # 共通(gene非依存)ゲノム変異ノードの三つ組。bulk の共有ファイル/自己完結モードで共用。
  def self.genomic_variant_triples(vid, v)
    cn = v['chr'].to_s.sub(/\Achr/, '')
    vtype = v['ref'].to_s.length == 1 && v['alt'].to_s.length == 1 ? 'gvo:SNV' : 'gvo:Variant'
    ["a #{vtype}",
     "dcterms:identifier #{Jogo.lit(vid)}",
     "gvo:ref #{Jogo.lit(v['ref'])}",
     "gvo:alt #{Jogo.lit(v['alt'])}",
     "faldo:location [ a faldo:ExactPosition ; faldo:position #{v['pos'].to_i} ; " \
     "faldo:reference hco:#{cn}\\/GRCh38 ]"]
  end

  attr_reader :gene_key, :mane_status, :genomic

  # siblings    : 同一遺伝子の別 region の gene URI (CURIE) 配列。rdfs:seeAlso に。
  # emit_common : true なら共通ゲノム変異ノードを当ファイルにインライン(自己完結=テスト用)。
  #               false なら出力せず @genomic に収集 (bulk で共有ファイルへ一度だけ出す)。
  def initialize(data, siblings: [], emit_common: true)
    @d = data
    @m = data.fetch('maneinfo')
    @symbol = @m['symbol']
    @region = @m['regionname5000']
    @mane_status = @m['mane_status']
    # gene_key: MANE Select は symbol、それ以外(MANE Plus Clinical)は regionname。
    # gene/haplotype/variant の全 URI をこのキーでスコープし region 間衝突を防ぐ。
    @gene_key = @mane_status == 'MANE Select' ? @symbol : @region
    @siblings = siblings || []
    @emit_common = emit_common
    @genomic = {}
    @out = +''
  end

  def to_ttl
    build_variant_index
    @genomic = {}
    @out = +''
    @out << header
    @out << "# Generated by jogo2ttl.rb — gene #{@symbol} (#{@region})\n\n"
    build_gene
    build_haplotypes
    build_variants
    @out
  end

  private

  # --- リテラル/URI ヘルパ ---------------------------------------------------
  def str(s) = Jogo.lit(s)

  def bool(v) = v.to_i.zero? ? 'false' : 'true'

  def present?(x) = !x.nil? && x.to_s != ''

  # ゼロパディング除去: a0001→a1, a0000→a0, c0010→c10, a0001c0001→a1c1
  def unpad(id) = id.to_s.gsub(/([A-Za-z])(\d+)/) { "#{Regexp.last_match(1)}#{Regexp.last_match(2).to_i}" }

  def hap_uri(padded) = "jogo_hap:#{@gene_key}:#{unpad(padded)}"

  def var_uri(vid) = "jogo_var:#{@gene_key}:#{vid}"

  # 複合id ("a0005c0005t0066g0097") から指定レベルのトークン ("t0066") を取り出す
  def extract_token(composite, level)
    composite.to_s.scan(/([actg])(\d+)/) do |l, n|
      return "#{l}#{n}" if l == level
    end
    nil
  end

  def chrnum(chr) = chr.to_s.sub(/\Achr/, '')

  def var_id(v) = "#{chrnum(v['chr'])}-#{v['pos']}-#{v['ref']}-#{v['alt']}"

  def emit(subject, triples)
    @out << self.class.block(subject, triples)
  end

  def header = self.class.prefix_header

  # --- variant インデックス (hap → variant のリンクを先に構築) ---------------
  def build_variant_index
    @variants = {}
    @hap_vars = Hash.new { |h, k| h[k] = [] }
    VARIANT_TABLES.each do |sec, level, col|
      (@d[sec] || []).each do |v|
        vid = var_id(v)
        @variants[vid] ||= v
        (v[col] || '').split(',').reject(&:empty?).each do |hid|
          tok = extract_token(hid, level)
          @hap_vars[hap_uri(tok)] << vid if tok
        end
      end
    end
  end

  # --- 遺伝子領域 -----------------------------------------------------------
  def build_gene
    t = []
    t << 'a jogo:GeneRegion'
    t << "rdfs:label #{str(@symbol)}"
    t << "dcterms:description #{str(@m['name'])}" if present?(@m['name'])
    t << "jogo:regionName #{str(@region)}"
    t << "jogo:maneStatus #{str(@m['mane_status'])}" if present?(@m['mane_status'])
    t << "jogo:maneVersion #{str(@m['ver'])}" if present?(@m['ver'])
    t << "foaf:page <#{BASE}/gene?genename=#{@region}>"
    GENE_HAP_PRED.each do |sec, pred|
      rows = @d[sec] || []
      next if rows.empty?

      idcol = SUMMARY.find { |x| x[0] == sec }[2]
      uris = rows.map { |r| hap_uri(r[idcol]) }.uniq
      t << "jogo:#{pred} #{uris.join(', ')}"
    end
    t << "rdfs:seeAlso hgnc:#{@m['hgncid']}" if present?(@m['hgncid'])
    ens = @m['ensemblid'].to_s.sub(/\.\d+\z/, '')
    t << "rdfs:seeAlso ensembl:#{ens}" unless ens.empty?
    @siblings.each { |sib| t << "rdfs:seeAlso #{sib}" }  # 同一遺伝子の別 region
    cn = chrnum(@m['chr'])
    t << "faldo:location [ a faldo:Region ;\n        " \
         "faldo:begin [ a faldo:ExactPosition ; faldo:position #{@m['start'].to_i} ; " \
         "faldo:reference hco:#{cn}\\/GRCh38 ] ;\n        " \
         "faldo:end [ a faldo:ExactPosition ; faldo:position #{@m['end'].to_i} ; " \
         "faldo:reference hco:#{cn}\\/GRCh38 ] ]"
    emit("jogo_gene:#{@gene_key}", t)
  end

  # --- ハプロタイプ (個別 A/C/T/G + 複合 AC/ACT/ACTG) -----------------------
  def build_haplotypes
    SUMMARY.each do |sec, cls, idcol, lencol, indiv|
      (@d[sec] || []).each do |r|
        padded = r[idcol]
        s = hap_uri(padded)
        t = []
        t << "a jogo:#{cls}"
        t << "dcterms:identifier #{str("#{@gene_key}:#{unpad(padded)}")}"
        t << "skos:altLabel #{str("#{@gene_key}:#{padded}")}"
        t << "foaf:page <#{BASE}/haplotype_detail?hapid=#{padded}&regionname=#{@region}>"
        t << "jogo:sequenceLength #{r[lencol].to_i}" if lencol && r[lencol]
        indiv.each { |pred, col| t << "jogo:#{pred} #{hap_uri(r[col])}" if present?(r[col]) }
        t << "jogo:totalCount #{r['total'].to_i}" if r.key?('total')
        t << "jogo:isGRCh38Reference #{bool(r['grch38p0'])}" if r.key?('grch38p0')
        t << "jogo:isCHM13Reference #{bool(r['chm13v2'])}" if r.key?('chm13v2')
        t << "jogo:isMANEReference #{bool(r['MANEref'])}" if r.key?('MANEref')
        vs = (@hap_vars[s] || []).uniq
        t << "jogo:hasVariant #{vs.map { |vid| var_uri(vid) }.join(', ')}" unless vs.empty?
        pop_freqs(r).each do |pop, cnt|
          t << "jogo:frequency [ a jogo:PopulationFrequency ; " \
               "jogo:population jogo_pop:#{pop} ; jogo:alleleCount #{cnt} ]"
        end
        emit(s, t)
      end
    end
  end

  # *_total 列 (5 super-population と 全体 total を除く) を集団頻度に
  def pop_freqs(r)
    r.keys.select { |k| k.end_with?('_total') }
     .map { |k| k.sub('_total', '') }
     .reject { |p| SUPERPOPS.include?(p) }
     .map { |p| [p, r["#{p}_total"].to_i] }
     .select { |_, c| c.positive? }
  end

  # --- variant (a/c/t/g variants を var_id で重複排除) -----------------------
  def build_variants
    @variants.each do |vid, v|
      # gene固有ノード: 頻度・snpEff注釈・hgvs。共通ノードへ jogo:genomicVariant。
      gene = ['a jogo:Variant']
      gene << "jogo:alleleCount #{v['ac'].to_i}" if v.key?('ac')
      gene << "jogo:alleleNumber #{v['an'].to_i}" if v.key?('an')
      gene << "jogo:alleleFrequency #{v['af']}" if present?(v['af'])
      gene << "jogo:annotation #{str(v['snpeff_annotation'])}" if present?(v['snpeff_annotation'])
      gene << "jogo:annotationImpact #{str(v['snpeff_annotation_impact'])}" if present?(v['snpeff_annotation_impact'])
      gene << "jogo:hgvsC #{str(v['snpeff_hgvs_c'])}" if present?(v['snpeff_hgvs_c'])
      gene << "jogo:hgvsP #{str(v['snpeff_hgvs_p'])}" if present?(v['snpeff_hgvs_p'])
      gene << "jogo:genomicVariant jogo_var:#{vid}"  # 順向き(遺伝子別→共通)
      emit(var_uri(vid), gene)

      # 逆向き(共通→遺伝子別)リンク。主語は共通ノードだが gene 別なので重複しない。
      emit("jogo_var:#{vid}", ["jogo:hasGeneVariant #{var_uri(vid)}"])

      # 共通(gene非依存)ノード: emit_common=true ならインライン、false なら @genomic に収集。
      if @emit_common
        emit("jogo_var:#{vid}", self.class.genomic_variant_triples(vid, v))
      else
        @genomic[vid] = v
      end
    end
  end
end

# === main (単一遺伝子テスト用。本番は bin/regions2ttl.rb) =====================
if $PROGRAM_NAME == __FILE__
  opts = {}
  OptionParser.new do |o|
    o.banner = 'Usage: ruby bin/jogo2ttl.rb (--genename NAME | --regionname NAME) [--out PATH]'
    o.on('--genename NAME') { |x| opts[:genename] = x }
    o.on('--regionname NAME') { |x| opts[:regionname] = x }
    o.on('--out PATH') { |x| opts[:out] = x }
  end.parse!(ARGV)

  data =
    if opts[:genename] then Jogo.fetch_by_genename(opts[:genename])
    elsif opts[:regionname] then Jogo.fetch_by_regionname(opts[:regionname])
    else abort('specify --genename or --regionname')
    end

  builder = TtlBuilder.new(data) # emit_common: true で自己完結ファイル
  ttl = builder.to_ttl
  if opts[:out] == '-'
    $stdout.write(ttl)
  else
    require 'fileutils'
    # 本番と同様 output/<mane_status>/<gene_key>.ttl に振り分け (ディレクトリは自動作成)
    status_dir = builder.mane_status.to_s.gsub(/\s+/, '_')
    out = opts[:out] || File.expand_path("../output/#{status_dir}/#{builder.gene_key}.ttl", __dir__)
    FileUtils.mkdir_p(File.dirname(out))
    File.write(out, ttl)
    warn "wrote #{out} (#{ttl.bytesize} bytes)"
  end
end
