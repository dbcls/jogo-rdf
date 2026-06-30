#!/usr/bin/env ruby
# frozen_string_literal: true
#
# regions2ttl.rb — JoGo の全 region をループし、遺伝子(region)毎の ttl を出力する
#
#   region list : https://jogo.csml.org/api/v1/regions/   (cursor ページング)
#   region data : https://jogo.csml.org/genicregion?regionname=...&format=json
#
# 仕様:
#   - gene_key = MANE Select は symbol / MANE Plus Clinical は regionname。
#     gene/haplotype/variant の全 URI をこのキーでスコープ (jogo2ttl.rb 参照)。
#   - rdfs:seeAlso は MANE Plus Clinical の ttl にのみ出力し、対応する MANE Select
#     の gene を指す (Select↔Plus Clinical 間のみ。Plus同士は張らない)。
#     Select の ttl には seeAlso を出さない (Select 単独 graph で完結する運用のため)。
#   - 変異の重複排除: gene固有ノード(頻度/注釈)は各遺伝子 ttl に、共通(gene非依存)
#     ゲノム変異ノード(class/ref/alt/location)は output/genomic_variants.ttl に
#     全遺伝子で「1回だけ」出力 (グローバル dedup)。両者は jogo:genomicVariant /
#     jogo:hasGeneVariant で双方向リンク。
#   - 出力レイアウト:
#       output/MANE_Select/<symbol>.ttl
#       output/MANE_Plus_Clinical/<regionname>.ttl
#       output/genomic_variants.ttl                 (共通ゲノム変異; 全遺伝子で共有)
#
# 使い方 (プロジェクト直下から):
#   ruby bin/regions2ttl.rb                  # 全件 (19287 region; 本番)
#   ruby bin/regions2ttl.rb --gene BRAF      # 1遺伝子(全region)だけ (テスト)
#   ruby bin/regions2ttl.rb --limit 20       # 先頭20 region で試走
#   オプション: --out-dir DIR, --sleep 0.2 (取得間ウェイト), --resume (既存ttlはスキップ)

require_relative 'jogo2ttl'
require 'fileutils'
require 'optparse'
require 'set'

REGIONS_API = 'https://jogo.csml.org/api/v1/regions/'

module Regions
  module_function

  def fetch_all(limit: 1000)
    all = []
    cursor = 0
    loop do
      res = Net::HTTP.get_response(URI("#{REGIONS_API}?limit=#{limit}&cursor=#{cursor}"))
      raise "HTTP #{res.code} for regions API" unless res.is_a?(Net::HTTPSuccess)

      body = JSON.parse(res.body)
      rows = body['data'] || []
      all.concat(rows)
      nc = body.dig('meta', 'pagination', 'next_cursor')
      break if rows.empty? || rows.size < limit || nc.nil? || nc == cursor

      cursor = nc
    end
    all
  end

  # region list の1行から gene_key を決める (jogo2ttl.rb の TtlBuilder と同一規則)
  def gene_key(region)
    region['mane_status'] == 'MANE Select' ? region['genename'] : region['regionname']
  end

  def status_dir(region)
    region['mane_status'].to_s.gsub(/\s+/, '_')
  end
end

# === main ====================================================================
if $PROGRAM_NAME == __FILE__
  opts = { out_dir: File.expand_path('../output', __dir__), sleep: 0.0 }
  OptionParser.new do |o|
    o.banner = 'Usage: ruby regions2ttl.rb [--gene SYMBOL] [--limit N] [--out-dir DIR] [--sleep SEC] [--resume]'
    o.on('--gene SYMBOL', 'この genename の region だけ処理') { |x| opts[:gene] = x }
    o.on('--limit N', Integer, '先頭 N region だけ処理') { |x| opts[:limit] = x }
    o.on('--out-dir DIR', '出力ルート (既定 プロジェクト直下 output/)') { |x| opts[:out_dir] = x }
    o.on('--sleep SEC', Float, '各 region 取得間のウェイト秒') { |x| opts[:sleep] = x }
    o.on('--resume', '既存 ttl が在れば取得・生成をスキップ') { opts[:resume] = true }
  end.parse!(ARGV)

  warn 'fetching region list ...'
  regions = Regions.fetch_all
  warn "#{regions.size} regions, #{regions.map { |r| r['genename'] }.uniq.size} genes"

  by_gene = regions.group_by { |r| r['genename'] }

  targets = regions
  targets = targets.select { |r| r['genename'] == opts[:gene] } if opts[:gene]
  targets = targets.first(opts[:limit]) if opts[:limit]
  warn "target regions: #{targets.size}"

  # 共通(gene非依存)ゲノム変異ノードは全遺伝子で1回だけ output/genomic_variants.ttl に出す。
  # seen で vid をグローバル dedup。--resume 時は既存ファイルから seen を復元して追記。
  shared_path = File.join(opts[:out_dir], 'genomic_variants.ttl')
  seen = Set.new
  if opts[:resume] && File.exist?(shared_path)
    File.foreach(shared_path) { |l| (m = l.match(/^jogo_var:(\S+)\s*$/)) && (seen << m[1]) }
    shared = File.open(shared_path, 'a')
    warn "resume: #{seen.size} genomic variants already in #{shared_path}"
  else
    FileUtils.mkdir_p(opts[:out_dir])
    shared = File.open(shared_path, 'w')
    shared << TtlBuilder.prefix_header
    shared << "# 共通(gene非依存)ゲノム変異ノード。各遺伝子 ttl の jogo:Variant から jogo:genomicVariant で参照。\n\n"
  end

  ok = 0
  err = 0
  skip = 0
  targets.each_with_index do |r, i|
    key = Regions.gene_key(r)
    dir = File.join(opts[:out_dir], Regions.status_dir(r))
    path = File.join(dir, "#{key}.ttl")
    if opts[:resume] && File.exist?(path)
      skip += 1
      next # 既存 gene ttl の変異は前回 shared に出済み
    end

    begin
      # rdfs:seeAlso: MANE Plus Clinical の時だけ、対応する MANE Select の gene を指す。
      # (Select の ttl には出さない / Plus 同士も張らない)
      siblings =
        if r['mane_status'] == 'MANE Select'
          []
        else
          sel = by_gene[r['genename']].find { |s| s['mane_status'] == 'MANE Select' }
          sel ? ["jogo_gene:#{Regions.gene_key(sel)}"] : []
        end
      data = Jogo.fetch_by_regionname(r['regionname'])
      builder = TtlBuilder.new(data, siblings: siblings, emit_common: false)
      FileUtils.mkdir_p(dir)
      File.write(path, builder.to_ttl)
      # 未出の共通ゲノム変異だけ shared に追記 (グローバル dedup)
      builder.genomic.each do |vid, v|
        next unless seen.add?(vid)

        shared << TtlBuilder.block("jogo_var:#{vid}", TtlBuilder.genomic_variant_triples(vid, v))
      end
      ok += 1
      warn "[#{i + 1}/#{targets.size}] #{r['mane_status']}  #{path}"
    rescue StandardError => e
      err += 1
      warn "[#{i + 1}/#{targets.size}] ERROR #{r['regionname']}: #{e.class} #{e.message}"
    end
    sleep(opts[:sleep]) if opts[:sleep].positive?
  end

  shared.close
  warn "done: ok=#{ok} skip=#{skip} err=#{err}  genomic variants=#{seen.size}"
  warn "  -> gene ttl: #{opts[:out_dir]}/<status>/  +  shared: #{shared_path}"
end
