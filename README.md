# jogo-rdf

## コンバート

```sh
# 1) 全遺伝子・全 region → 遺伝子別 ttl ＋ 共通ゲノム変異 ttl (変異は重複排除)
ruby bin/regions2ttl.rb --sleep 0.2 --resume     # 19,287 region。--resume で中断再開可

# 2) sample → output/samples.ttl
ruby bin/samples2ttl.rb
```

## 出力
* MANE plus clinicalはURIやファイル名がgene reagion name

```
output/
├── MANE_Select/<symbol>.ttl 
├── MANE_Plus_Clinical/<regionname>.ttl
├── genomic_variants.ttl                     # gvo, ref, altなどの共通部分
└── samples.ttl                              # sample, population
```

## 動作確認 (単一遺伝子・テスト用)

```sh
ruby bin/jogo2ttl.rb --genename ALDH2                                 # → output/MANE_Select/ALDH2.ttl
ruby bin/jogo2ttl.rb --regionname BRAF_chr7_140714337_140929929       # → output/MANE_Plus_Clinical/…
ruby bin/jogo2ttl.rb --genename ALDH2 --out -                         # 標準出力
ruby bin/regions2ttl.rb --gene BRAF                                   # 1遺伝子の全 region (本番と同じ分割出力)
```

## スクリプト一覧

| スクリプト | 役割 | 出力 |
|---|---|---|
| `bin/regions2ttl.rb` | **本番**: 全 region 一括変換（変異 dedup） | `output/<status>/*.ttl` ＋ `output/genomic_variants.ttl` |
| `bin/samples2ttl.rb` | 集団 + sample 変換（`templates/population.ttl` ＋ sample） | `output/samples.ttl` |
| `bin/jogo2ttl.rb` | 単一遺伝子の確認用（`TtlBuilder` 本体） | `output/<status>/<gene_key>.ttl`（自己完結） |

## 依存

- Ruby 3.x（標準ライブラリのみ）
- ネットワーク（全データを API から取得）
