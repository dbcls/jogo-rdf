# JoGo RDF コンバータ

JoGo Haplotype Database の API を、rdf-config の `model.yaml` に準拠した RDF/Turtle に
変換する。本番は **全遺伝子・全 region を一括変換** ＋ **sample を変換** の2ステップ。
Ruby 標準ライブラリのみ・依存 gem なし。入力は常に API（ローカル JSON は使わない）。

## 本番フロー (プロジェクト直下から)

```sh
# 1) 全遺伝子・全 region → 遺伝子別 ttl ＋ 共通ゲノム変異 ttl (変異は重複排除)
ruby bin/regions2ttl.rb --sleep 0.2 --resume     # 19,287 region。--resume で中断再開可

# 2) sample → output/samples.ttl
ruby bin/samples2ttl.rb
```

## 出力レイアウト

```
output/
├── MANE_Select/<symbol>.ttl                # 例 ALDH2.ttl   (遺伝子別: gene/haplotype/頻度/gene固有variant)
├── MANE_Plus_Clinical/<regionname>.ttl     # 例 BRAF_chr7_….ttl
├── genomic_variants.ttl                     # 共通(gene非依存)ゲノム変異。全遺伝子で1ファイルに集約
└── samples.ttl                              # 集団/スーパー集団(1000G) + sample (単体で完結)
```

`mane_status` 別にディレクトリを分割（エンドポイント/グラフを status で分ける運用を想定）。
ディレクトリは自動作成される。`output/samples.ttl` は読み込み用テンプレート
`templates/population.ttl`（集団定義）に sample を足した**完全版**で、単体でロードできる。

## 重複排除と variant の設計

variant を2ノードに分け、座標ベースの共通部分を1ファイルに集約して重複を避ける。

- **共通ノード** `jogo_var:<chr-pos-ref-alt>` … `a gvo:SNV` / `gvo:ref` / `gvo:alt` /
  `faldo:location`。座標ベースで全遺伝子共通 → `output/genomic_variants.ttl` に **1回だけ** 出力。
- **遺伝子別ノード** `jogo_var:<gene_key>:<chr-pos-ref-alt>` … `a jogo:Variant` /
  頻度(`alleleCount/Number/Frequency`) / snpEff注釈 / hgvs。遺伝子毎に異なるので各遺伝子 ttl に。
- **双方向リンク**: 遺伝子別 `jogo:genomicVariant →` 共通 / 共通 `jogo:hasGeneVariant →` 遺伝子別。
  同一物理変異を複数遺伝子が参照しても、共通ノードがハブとして各遺伝子レコードを束ねる。

## URI スコープ (gene_key)

同一遺伝子が複数 region(=複数 MANE 転写産物)を持つ34遺伝子でも衝突しないよう、
全 URI を **gene_key** でスコープする。`gene_key = MANE Select→symbol / MANE Plus Clinical→regionname`。

- gene `jogo_gene:<gene_key>` / haplotype `jogo_hap:<gene_key>:a1` / 遺伝子別variant `jogo_var:<gene_key>:…`
- haplotype の `dcterms:identifier` はゼロパディング除去(`a0001→a1`)、`skos:altLabel` はパディング付き。階層部分はコロン無し(`a1c1t1g4`)。
- `rdfs:seeAlso`: MANE Plus Clinical の ttl にのみ、対応する MANE Select gene を指す（Select 側には出さない）。

## 動作確認 (単一遺伝子・テスト用)

`bin/jogo2ttl.rb` は1遺伝子だけを **自己完結 ttl**（共通ノードもインライン）で出力する確認用ツール。

```sh
ruby bin/jogo2ttl.rb --genename ALDH2                                 # → output/MANE_Select/ALDH2.ttl
ruby bin/jogo2ttl.rb --regionname BRAF_chr7_140714337_140929929       # → output/MANE_Plus_Clinical/…
ruby bin/jogo2ttl.rb --genename ALDH2 --out -                         # 標準出力
ruby bin/regions2ttl.rb --gene BRAF                                   # 1遺伝子の全 region (本番と同じ分割出力)
```

出力は本番同様 `output/<mane_status>/<gene_key>.ttl` に振り分けられる（自動 mkdir）。

## sample の絞り込み

`bin/samples2ttl.rb` は `groupnames` に `JoGo_JPT_HPRC` / `JoGo_JPT` を含む sample のみを対象
（`KEEP_GROUPS` で指定）。`templates/population.ttl` に sample を足した `output/samples.ttl`
を毎回まるごと再生成する。

## 関連する静的ファイル

読み込み用テンプレート (`templates/`):

- `population.ttl` — 集団/スーパー集団(1000G準拠)の **curated マスタ**。手で編集する元データで、`samples2ttl.rb` が sample を足して `output/samples.ttl` を生成する。

rdf-config 設定・モデル (`sample/`):

- `ontology.ttl` — JoGo オントロジー（クラス/プロパティ定義）
- `prefix.yaml` — prefix 定義
- `model.yaml` — rdf-config モデル

ロード時は `output/**.ttl`（gene / genomic_variants / samples=集団含む）＋ `sample/ontology.ttl` を併せて使う（`templates/population.ttl` は samples.ttl に含まれるため別途ロード不要）。

## スクリプト一覧

| スクリプト | 役割 | 出力 |
|---|---|---|
| `bin/regions2ttl.rb` | **本番**: 全 region 一括変換（変異 dedup） | `output/<status>/*.ttl` ＋ `output/genomic_variants.ttl` |
| `bin/samples2ttl.rb` | 集団 + sample 変換（`templates/population.ttl` ＋ sample） | `output/samples.ttl` |
| `bin/jogo2ttl.rb` | 単一遺伝子の確認用（`TtlBuilder` 本体） | `output/<status>/<gene_key>.ttl`（自己完結） |

## 依存

- Ruby 3.x（標準ライブラリのみ）
- ネットワーク（全データを API から取得）
