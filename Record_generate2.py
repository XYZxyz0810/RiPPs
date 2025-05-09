from collections import defaultdict
import os
import logging
from Bio import SeqIO


def process_prodigal_faa_output(faa_path, output_dir, pre_min, pre_max):
    """
    处理 Prodigal 预测结果，生成两个输出文件：
    1. output_file: 记录所有 ORF 信息
    2. filtered_output_file: 仅记录符合长度阈值的 ORF

    :param faa_path: Prodigal 生成的 FAA 文件路径
    :param gbk_path: Prodigal 生成的 GBK 文件路径
    :param output_dir: 输出目录
    :param pre_min: 过滤 ORF 的最小长度
    :param pre_max: 过滤 ORF 的最大长度
    :return: 生成的 output_file 和 filtered_output_file 路径
    """
    CDS_file = os.path.join(output_dir, f"Temp_CDS_info.txt")
    peptide_file = os.path.join(output_dir, f"Temp_Peptide_info.txt")

    total_sequences = 0
    filtered_count = 0

    with open(faa_path, "r") as faa_f, \
            open(CDS_file, "w") as CDS_f, \
            open(peptide_file, "w") as peptide_f:

        # 写入表头
        CDS_f.write("accession\tstart\tend\tdirection\tlength\n")
        peptide_f.write("accession\tstart\tend\tdirection\tlength\tsequence\n")

        for record in SeqIO.parse(faa_f, "fasta"):
            total_sequences += 1
            try:
                # 解析 FASTA 头部信息
                parts = record.description.split(" # ")
                if len(parts) < 4:
                    logging.warning(f"Skipping malformed header: {record.description}")
                    continue

                accession = parts[0].lstrip(">")  # Accession ID
                start, end, direction = parts[1], parts[2], parts[3]
                sequence = str(record.seq)
                sequence_length = len(sequence)

                # 写入完整输出
                CDS_f.write(f"{accession}\t{start}\t{end}\t{direction}\t{sequence_length}\n")

                # 过滤符合长度要求的 ORF
                if pre_min <= sequence_length <= pre_max:
                    peptide_f.write(
                        f"{accession}\t{start}\t{end}\t{direction}\t{sequence_length}\t{sequence}\n")
                    filtered_count += 1

            except Exception as e:
                logging.error(f"Error processing sequence {record.id}: {e}")

    logging.info(f"Total ORFs: {total_sequences}, Filtered: {filtered_count}")
    print(f"Processed Prodigal output saved to {CDS_file}")
    print(f"Filtered sequences saved to {peptide_file}")

    return CDS_file, peptide_file


def comb_CDS_and_hmmscan_files(out_dir, CDS_file, parsed_hmmscan_results):
    """
    处理当前目录下 out_dir 中的 1_processed_faa_info.txt 和 parsed_hmmscan_results.txt 文件，
    以 accession 为 id 建立字典，并将相关信息组合输出成 combined_results.txt 文件。
    如果文件已存在，则覆盖。

    :param out_dir: 输出目录
    :param CDS_file: Prodigal 预测的 faa 文件路径
    :param hmmscan_file: Hmmscan 结果文件路径
    """
    # 设置文件路径
    CDS_prediction_file = os.path.join(out_dir, "Temp_combined_info.txt")

    # 建立 hmmscan_results 字典 (以 accession 为键，值为匹配的所有结果)
    hmmscan_dict = {}
    with open(parsed_hmmscan_results, "r") as f:
        for line in f:
            # 假设 parsed_hmmscan_results.txt 格式: accession\tPfam/Hmm\tname\tdescription\tE-value
            fields = line.strip().split("\t")
            if len(fields) >= 5:
                accession = fields[0]
                pfam_accession = fields[1]
                hit_name = fields[2]
                hit_description = fields[3]
                evalue = fields[4]

                # 如果 accession 已经存在于字典中，则追加结果
                if accession not in hmmscan_dict:
                    hmmscan_dict[accession] = []
                hmmscan_dict[accession].append([pfam_accession, hit_name, hit_description, evalue])

    # 打开 hmmscan_results.txt 文件，使用 'w' 模式覆盖文件内容并写入表头
    with open(CDS_prediction_file, "w") as combined_out:
        # 解析 faa_info 文件并组合信息
        with open(CDS_file, "r") as f:
            for line in f:
                # 假设格式: accession start end direction length
                fields = line.strip().split("\t")
                if len(fields) >= 5:
                    accession = fields[0]
                    start = fields[1]
                    end = fields[2]
                    direction = fields[3]
                    length = fields[4]

                    # 如果 hmmscan_dict 中存在此 accession，则遍历所有相关的结果
                    if accession in hmmscan_dict:
                        for result in hmmscan_dict[accession]:
                            pfam_accession, hit_name, hit_description, evalue = result
                            # 输出到 combined_results.txt
                            combined_out.write(
                                f"{accession}\t{start}\t{end}\t{direction}\t{length}\t{pfam_accession}\t{hit_name}\t{hit_description}\t{evalue}\n")
                    else:
                        # 如果没有找到相关信息，则输出 'none'
                        combined_out.write(
                            f"{accession}\t{start}\t{end}\t{direction}\t{length}\tnone\tnone\tnone\tnone\n")

    print(f"Combined results have been saved to {CDS_prediction_file}")

    return CDS_prediction_file


def generate_cluster_prediction_tables(out_dir, CDS_prediction_file):
    """
    根据 accession 中 _ 前的内容对数据进行分组，并生成多个表格，
    对相同 accession 的行，将重复行中的 accession, start, end, direction, length 列替换为空字符。

    :param out_dir: 输出目录
    :param combined_output_file: combined_results.txt 文件路径
    """
    # 使用 defaultdict 来存储按 accession 前缀分组的结果
    grouped_data = defaultdict(list)
    seen_accessions = set()  # 记录已处理的 accession
    # 读取 combined_results.txt 文件并跳过第一行
    with open(CDS_prediction_file, "r") as f:
        first_line = True
        for line in f:
            if first_line:
                first_line = False
                continue  # 跳过第一行（表头）

            fields = line.strip().split("\t")
            if len(fields) >= 9:
                accession = fields[0]
                start = fields[1]
                end = fields[2]
                direction = fields[3]
                length = fields[4]
                pfam_accession = fields[5]
                hit_name = fields[6]
                hit_description = fields[7]
                evalue = fields[8]

                # 提取 accession 中 _ 前的内容作为分组依据
                prefix = "_".join(accession.split("_")[:-1])

                # 检查是否是该 accession 的首次出现
                if accession not in seen_accessions:
                    seen_accessions.add(accession)
                    grouped_data[prefix].append(
                        [accession, start, end, direction, length, pfam_accession, hit_name, hit_description, evalue])
                else:
                    # 对重复出现的 accession 行，使用空字符替代某些列
                    grouped_data[prefix].append(
                        ["", "", "", "", "", pfam_accession, hit_name, hit_description, evalue]
                    )

    # 为每个 prefix 生成一个对应的文件并返回文件路径
    Cluster_prediction_table = os.path.join(out_dir, "Final_BGC_table.txt")
    with open(Cluster_prediction_table, "w") as out:
        for prefix, data in grouped_data.items():
            # 写入分隔线和表头
            out.write("\n" + "=" * 50 + f"\nTable for prefix: {prefix}\n" + "=" * 50 + "\n")

            # 定义列宽（根据最长的内容自动调整）
            col_widths = [20, 10, 10, 10, 10, 15, 20, 30, 10]  # 每列的固定宽度

            # 格式化表头
            header = ["accession", "start", "end", "direction", "length", "Pfam/Hmm", "name", "description", "E-value"]
            out.write(" ".join(f"{col:<{col_widths[i]}}" for i, col in enumerate(header)) + "\n")

            # 写入数据
            for row in data:
                out.write(" ".join(f"{field:<{col_widths[i]}}" for i, field in enumerate(row)) + "\n")

        print(f"Grouped results have been saved to {Cluster_prediction_table}")

    # 返回文件路径
    return Cluster_prediction_table

def generate_precursor_fasta(out_dir, peptide_file):
    """
    处理 peptide 文件，读取文件内容，并根据 accession 为索引建立字典，
    然后以指定格式写入输出文件 Peptide_file.txt。

    :param out_dir: 输出目录
    :param peptide_file: peptide 文件路径
    """
    peptide_dict = {}

    # 读取 peptide_file 内容并根据 accession 建立字典
    with open(peptide_file, "r") as f:
        header = f.readline()  # 跳过表头
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 6:
                accession = fields[0]
                start = fields[1]
                end = fields[2]
                direction = fields[3]
                length = fields[4]
                sequence = fields[5]

                # 去除序列末尾的 '*' 和换行符
                sequence = sequence.rstrip("*").replace("\n", "")

                # 将 accession 作为键，相关数据作为值存储在字典中
                peptide_dict[accession] = {
                    "start": start,
                    "end": end,
                    "direction": direction,
                    "length": length,
                    "sequence": sequence
                }

    # 设置输出文件路径
    peptide_fasta_file = os.path.join(out_dir, "Final_peptide_candidate.txt")

    # 根据 accession 为索引的顺序将数据写入输出文件
    with open(peptide_fasta_file, "w") as out:

        # 写入 peptide 数据
        for accession, data in peptide_dict.items():
            start = data["start"]
            end = data["end"]
            direction = data["direction"]
            length = data["length"]
            sequence = data["sequence"]

            out.write(f">{accession}\t{start}\t{end}\t{direction}\t{length}\n{sequence}\n")

    print(f"Peptide data has been saved to {peptide_fasta_file}")

    return peptide_fasta_file


if __name__ == "__main__":
    # 获取当前目录下的 out_dir
    current_dir = os.getcwd()
    out_dir = os.path.join(current_dir, "out_dir")

    # 假设 output_file 和 hmmscan_output 文件路径已定义
    #output_file = os.path.join(out_dir, "1_10kb_nucl_CDS_info.txt")
    #hmmscan_output = os.path.join(out_dir, "parsed_hmmscan_results.txt")

    # 运行文件处理和生成分组表
    #combined_output_file = comb_files(out_dir, output_file, hmmscan_output)
    #result_file_path = generate_tables(out_dir, combined_output_file)

    # 打开并读取 grouped_results.txt 文件内容
    #with open(result_file_path, "r") as result_file:
        #print(result_file.read())  # 读取并打印文件内容

    peptide_file = os.path.join(out_dir, "Peptide_info.txt")
    peptide_fasta_file = generate_precursor_fasta(out_dir, peptide_file)