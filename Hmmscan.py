import os
import subprocess
from Bio import SearchIO


def run_hmmscan(faa_path, hmm_db_path, output_dir):
    """
    运行 HMMScan 对 .faa 文件进行扫描

    :param faa_path: 输入的 .faa 文件路径
    :param hmm_db_path: HMM 数据库路径
    :param output_dir: 输出目录
    :return: HMMScan 结果文件路径
    """
    # 确保输出目录存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 设置 hmmscan 输出文件路径
    hmmscan_output = os.path.join(output_dir, "hmmscan_results.txt")

    # 运行 hmmscan
    hmmscan_cmd = [
        "hmmscan",
        "-E", "1e-5",  # E 值阈值
        "--domtblout", hmmscan_output,  # 结果保存路径
        hmm_db_path,  # HMM 数据库
        faa_path  # 输入 .faa 文件
    ]

    try:
        subprocess.run(hmmscan_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running hmmscan: {e}")
        return None

    return hmmscan_output  # 返回 hmmscan 输出结果文件的路径


def process_hmmscan_output(hmmscan_output, output_dir):
    """
    解析 hmmscan 的输出文件，提取结构域信息并保存

    :param hmmscan_output: HMMScan 结果文件路径
    :param output_dir: 输出目录
    """
    if not os.path.exists(hmmscan_output):
        print(f"Warning: {hmmscan_output} does not exist.")
        return

    # 解析并保存文件路径
    parsed_hmmscan_results = os.path.join(output_dir, "hmmscan_info.txt")

    with open(parsed_hmmscan_results, "w") as out_f:
        out_f.write("accession\tPfam/Hmm\tname\tdescription\tE-value\n")

        # 解析 hmmscan 输出文件
        for result in SearchIO.parse(hmmscan_output, "hmmscan3-domtab"):
            query_name = result.id  # 查询序列 ID

            for hit in result.hits:
                pfam_accession = hit.accession.split('.')[0]  # 提取 PFam 号
                target_name = hit.id
                target_description = hit.description
                evalue = hit.evalue  # 选择 e-value

                out_f.write(f"{query_name}\t{pfam_accession}\t{target_name}\t{target_description}\t{evalue}\n")

    print(f"Parsed results saved to {parsed_hmmscan_results}")
    return parsed_hmmscan_results  # 返回解析后的结果文件路径
