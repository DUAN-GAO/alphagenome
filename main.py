import argparse
import myvariant
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import variant_scorers


# 固定 API Key
API_KEY = "AIzaSyD5Kht8QzCPkHeJ456_Tf_eBWirtKhmaRU"


def rsid_to_variant_info(rsid, genome_build="hg19"):
    """
    输入 rsID，返回:
    chromosome='chr22', position=36201698, reference_bases='A', alternate_bases='C'
    仅适用 SNV
    """
    mv = myvariant.MyVariantInfo()
    out = mv.getvariant(rsid, fields="dbsnp")
    if not out or "dbsnp" not in out:
        raise ValueError(f"{rsid} 未找到相关 dbSNP 信息")

    dbsnp = out["dbsnp"]
    vartype = dbsnp.get("vartype", "").lower()
    if vartype != "snv":
        raise ValueError(f"{rsid} 不是单核苷酸变异 (SNV)")

    ref_base = dbsnp.get("ref")
    alt_base = dbsnp.get("alt")
    chrom = dbsnp.get("chrom")
    pos = dbsnp.get(genome_build, {}).get("start")

    if None in [ref_base, alt_base, chrom, pos]:
        raise ValueError(f"{rsid} 坐标或碱基信息不完整")

    return f"chr{chrom}", pos, ref_base, alt_base


def main(rsid):
    # 初始化模型
    print('API reached...')
    dna_model = dna_client.create(API_KEY)
    
    # # 下载并处理 GTF 文件（只需第一次）
    # gtf = pd.read_feather(
    #     "https://storage.googleapis.com/alphagenome/reference/gencode/"
    #     "hg38/gencode.v46.annotation.gtf.gz.feather"
    # )
    # gtf_transcripts = gene_annotation.filter_protein_coding(gtf)
    # gtf_transcripts = gene_annotation.filter_to_longest_transcript(gtf_transcripts)
    # transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcripts)

    # rsID 转换成 variant
    chrom, pos, ref, alt = rsid_to_variant_info(rsid)
    print(f"[INFO] {rsid} -> {chrom}:{pos} {ref}>{alt}")

    variant = genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt,
    )

    # 定义预测区间
    sequence_length = 2048
    interval = variant.reference_interval.resize(sequence_length)

    # 定义 scorer（CenterMaskScorer）
    scorer = variant_scorers.CenterMaskScorer(
        width=None,
        aggregation_type=variant_scorers.AggregationType.DIFF_SUM_LOG2,
        requested_output=dna_client.OutputType.RNA_SEQ,
    )

    # 打分
    score_result = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=[scorer],
        organism=dna_client.Organism.HOMO_SAPIENS,
    )

    # 保存到以 rsID 命名的 txt 文件
    output_file = f"{rsid}.txt"
    with open(output_file, "w") as f:
        f.write(str(score_result[0].var))

    print(f"[RESULT] 已保存 score_result[0].var 到 {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Query rsID and score variant with AlphaGenome.")
    parser.add_argument("--rsid", type=str, required=True, help="dbSNP rsID, e.g. rs5934683")

    args = parser.parse_args()
    main(args.rsid)
