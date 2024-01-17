import pandas as pd
import sys
from gtfparse import read_gtf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import fisher_exact
import seaborn as sns

mpl.rcParams["pdf.fonttype"]=42


def get_age_localization(location, gtf, age, output):
    location = pd.read_csv(location, sep=",")
    gtf = read_gtf(gtf)
    age = pd.read_csv(age, sep="\t")
    protein = gtf[gtf["gene_biotype"]=="protein_coding"]
    protein2gene = protein.set_index("protein_id")["gene_id"].to_dict()
    gene =  gtf[gtf["feature"]=="gene"]
    gene2chromosome = gene.set_index("gene_id")["seqname"]
    gene2age = age.set_index("Gene")["Branch"].to_dict()
    localization = location[["Localizations", "Protein_ID"]]
    localization["Gene_id"] = localization["Protein_ID"].map(protein2gene)
    localization["Age"] = localization["Gene_id"].map(gene2age)
    localization["Chr"] = localization["Gene_id"].map(gene2chromosome)
    localization[["Gene_id", "Protein_ID", "Chr", "Age", "Localizations"]].to_csv(output, sep="\t", index=False)

def extract_number(s):
    letter_part = ''.join(filter(str.isalpha, s))
    number_part = int(''.join(filter(lambda x: x.isdigit() or x == '-', s)))
    return letter_part, number_part


if __name__ == "__main__":
    location = sys.argv[1]
    gtf = sys.argv[2]
    age = sys.argv[3]
    output = sys.argv[4]
    #get_age_localization(location, gtf, age, output)

    branches = ["br-2", "br-1", "br0", "br1", "br2", "br3", "br4", "br5", "br6"]
    groups = ["Oldest", "Older", "Old", "Young"]
    # 开始进行计算、分析、绘图
    age_localization = pd.read_csv(output, sep="\t")
    age_localization = age_localization.dropna()
    age_localization = age_localization[age_localization["Chr"].isin(["2L", "2R", "3L","3R","4","X","Y","mitochondrion_genome"])]
    brs = list(set(age_localization["Age"]))
    sorted_brs = sorted(brs, key=extract_number)
    compartments = set(age_localization["Localizations"].str.cat(sep="|").split("|"))
    compartments.remove("Plastid")
    compartments = list(compartments)
    compartments.sort()
    br_compartment_array = []
    br_ratio_list = []
    for br in sorted_brs:
        br_tmp = []
        br_gene_num=len(set(age_localization[age_localization["Age"]==br]["Gene_id"]))
        all_gene_num = len(set(age_localization["Gene_id"]))
        #print(br, br_gene_num, all_gene_num)
        br_ratio = round(br_gene_num/all_gene_num,6)
        br_ratio_list.append(br_ratio)
        for i in compartments:
            i_br_compartment = set(age_localization[(age_localization["Localizations"].str.contains(i)) & (age_localization["Age"]==br)]["Gene_id"])
            #print(f"{br}\t{i}\t{len(i_br_compartment)}")
            br_tmp.append(len(i_br_compartment))
        
        br_compartment_array.append(br_tmp)
    br_compartment_array = np.array(br_compartment_array)
    br_compartment_pd = pd.DataFrame(br_compartment_array, index=sorted_brs, columns=compartments)
    br_compartment_pd.to_csv("del", sep="\t")
    sum_br_compartment = list(br_compartment_pd.sum(axis=0))
    br_compartment_array_expect = []
    for j in br_ratio_list:
        tmp_expect = []
        for i in sum_br_compartment:
            tmp_expect.append(int(i*j))
        br_compartment_array_expect.append(tmp_expect)
    br_compartment_expect_pd = pd.DataFrame(br_compartment_array_expect,index=sorted_brs,columns=compartments)
    br_compartment_expect_pd.to_csv("del", sep="\t")
    ei = (br_compartment_pd-br_compartment_expect_pd)/br_compartment_expect_pd

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    colours = ["#050f2c", "#003666", "#00aeff", "#3369e7", "#8e43e7", "#b84592", "#ff4f81", "#ff6c5f", "#ffc168"]
    ei.plot(kind='bar', ax=axs[0], color=colours)
    axs[0].legend(labels=[])
    # 添加垂直于x轴的线
    for i in range(1, ei.shape[0]):
        axs[0].axvline(x=i - 0.5, color='grey', linestyle='--', linewidth=1)
    axs[0].set_ylabel('Excess percentage')
    axs[0].axhline(0, color='black', linestyle='-', linewidth=0.5)
    br_ratio_list=[]
    br_compartment_array=[]

    oldest = sorted_brs[0]; older=sorted_brs[1]; old=sorted_brs[2]; young = sorted_brs[3:]
    new_sorted_brs = [[oldest], [older], [old], young]
    for br in new_sorted_brs:
        br_tmp = []
        br_gene_num=len(set(age_localization[age_localization["Age"].isin(br)]["Gene_id"]))
        all_gene_num = len(set(age_localization["Gene_id"]))
        br_ratio = round(br_gene_num/all_gene_num,6)
        br_ratio_list.append(br_ratio)
        for i in compartments:
            i_br_compartment = set(age_localization[(age_localization["Localizations"].str.contains(i)) & (age_localization["Age"].isin(br))]["Gene_id"])
            br_tmp.append(len(i_br_compartment))

        br_compartment_array.append(br_tmp)
    br_compartment_array = np.array(br_compartment_array)
    br_compartment_pd = pd.DataFrame(br_compartment_array, index=["Oldest", "Older", "Old", "Young"], columns=compartments)
    #br_compartment_pd.to_csv("del", sep="\t")
    #print(br_compartment_pd)
    
    # 进行fisher精确检验
    p_values = np.empty_like(br_compartment_pd.values, dtype=float)
    for col_idx, col in enumerate(br_compartment_pd.columns):
        for row_idx in range(br_compartment_pd.shape[0]):
            contingency_table = np.array([[br_compartment_pd.iloc[row_idx, col_idx], br_compartment_pd.iloc[row_idx, :].sum() - br_compartment_pd.iloc[row_idx, col_idx]],
                                          [br_compartment_pd.iloc[:, col_idx].sum() - br_compartment_pd.iloc[row_idx, col_idx], br_compartment_pd.values.sum() - br_compartment_pd.iloc[row_idx, :].sum() - br_compartment_pd.iloc[:, col_idx].sum() + br_compartment_pd.iloc[row_idx, col_idx]]])

            _, p_values[row_idx, col_idx] = fisher_exact(contingency_table)
            print(col_idx, col, p_values[row_idx, col_idx])
            
    log_p_values = -np.log10(p_values)
    # 绘制热图
    plt.figure(figsize=(10, 8))
    sns.heatmap(p_values, annot=False, cmap="coolwarm", fmt=".3f", cbar_kws={'label': 'P-value'}, yticklabels=["Oldest", "Older", "Old", "Young"], xticklabels = br_compartment_pd.columns)
    plt.title("Fisher Exact Test P-values Heatmap")
    
    plt.tight_layout()
    plt.savefig("/data/chuand/dm_subcellular/figure/dm_heatmap_4group.pdf")

    '''
    # 按照行进行求和，相当于对每一列求和
    sum_br_compartment = list(br_compartment_pd.sum(axis=0))
    br_compartment_array_expect = []
    for j in br_ratio_list:
        # j表示的是每一个branch上的基因占左右基因的比例
        tmp_expect = []
        for i in sum_br_compartment:
            # i表示的i细胞器上在所有br上定位的总和
            tmp_expect.append(int(i*j)) # 如果j表示的是br-2上的比例，则表示br-2上期望有多少个基因出现
        br_compartment_array_expect.append(tmp_expect)
    br_compartment_expect_pd = pd.DataFrame(br_compartment_array_expect,index=["Oldest", "Older", "Old", "Young"],columns=compartments)
    print(br_ratio_list)
    br_compartment_expect_pd.to_csv("del", sep="\t")
    ei = (br_compartment_pd-br_compartment_expect_pd)/br_compartment_expect_pd
    for i in range(1, ei.shape[0]):
        axs[1].axvline(x=i - 0.5, color='grey', linestyle='--', linewidth=1)
    # br_compartment_pd.to_csv("del", sep="\t")
    br_compartment_expect_pd.to_csv("del", sep="\t")
    print(ei)
    ei.plot(kind='bar', ax=axs[1], color=colours)
    #plt.subplots_adjust(bottom=0.5)
    legend = axs[1].legend(bbox_to_anchor=(0.5, 0.99), loc='upper center',ncols=3, columnspacing=0.3, handletextpad=0.2, fontsize=8)
    legend.get_frame().set_linewidth(0) 
    legend.get_frame().set_facecolor('none')
    legend._legend_box.align = "left"
    axs[1].set_ylabel('Excess percentage')
    axs[1].axhline(0, color='black', linestyle='-', linewidth=0.5)
    plt.tight_layout()
    plt.savefig("/data/chuand/dm_subcellular/figure/dm_compartment1.pdf")
    '''














