import os
import subprocess
import logging

# 配置日志
logging.basicConfig(filename="prodigal_processing.log", level=logging.INFO,
                    format="%(asctime)s - %(levelname)s - %(message)s")


def run_prodigal(input_fasta, output_dir):
    """
    运行 Prodigal 进行 ORF 预测

    :param input_fasta: 输入的 FASTA 文件路径
    :param output_dir: 输出目录
    :return: 生成的 FAA 和 GBK 文件路径
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    prodigal_faa = os.path.join(output_dir, f"Temp_sequence.faa")
    prodigal_gbk = os.path.join(output_dir, f"Temp_genebank.gbk")

    # 运行 Prodigal
    prodigal_cmd = [
        "prodigal",
        "-i", input_fasta,
        "-a", prodigal_faa,
        "-o", prodigal_gbk,
        "-f", "gbk",
        "-p", "meta"  # 以 meta 模式运行
    ]

    logging.info(f"Running Prodigal on {input_fasta}...")
    subprocess.run(prodigal_cmd, check=True)
    logging.info(f"Prodigal output saved to {prodigal_faa} and {prodigal_gbk}")

    return prodigal_faa, prodigal_gbk

def main():
    """
    运行完整流程：
    1. 运行 Prodigal 预测 ORF
    2. 处理 Prodigal 结果并进行长度筛选
    """
    # 直接在这里定义参数
    pre_min = 100  # 过滤 ORF 的最小长度
    pre_max = 1000  # 过滤 ORF 的最大长度
    fna_in_path = "GT2_LanC_10kb_nucl.fasta"  # 输入的 FASTA 文件路径
    out_dir = "output_dir2"  # 输出目录路径

    # 运行 Prodigal
    prodigal_faa, prodigal_gbk = run_prodigal(fna_in_path, out_dir)

    # 处理 Prodigal 输出，并获取生成的文件路径
    #output_file, filtered_output_file = process_prodigal_output(prodigal_faa, prodigal_gbk, out_dir, pre_min, pre_max)

    # 输出生成的文件路径
    #print(f"Processed output file: {output_file}")
    #print(f"Filtered output file: {filtered_output_file}")

    return prodigal_faa, prodigal_gbk


if __name__ == "__main__":
    main()

