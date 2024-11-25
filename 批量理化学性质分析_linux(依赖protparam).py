#!python
# coding: utf-8
# usage: python script infile filename
import re

from Bio import SeqIO
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By

file_name = input('请输入文件名:')
out_file = 'result.txt'

# 更换为自己的chromedriver.exe所在的位置
chrome_options = Options()
chrome_options.add_argument('--headless')
chrome_options.add_argument('-no-sandbox')
service = Service(r"D:\chromedriver\chromedriver.exe")
expasy = webdriver.Chrome(service=service, options=chrome_options)
expasy.get("https://web.expasy.org/protparam/")


class expasy_cal():
    # get physical and chemical parameters for a given protein sequence file based on web
    # https://web.expasy.org/protparam/
    def inputseq(seq):
        """input the protein sequence"""
        while True:
            if expasy.find_element(By.XPATH, '/html/body/main/div/form/textarea').is_displayed():
                expasy.find_element(By.XPATH, '/html/body/main/div/form/input[2]').click()  # 获取新网页
                expasy.find_element(By.XPATH, '/html/body/main/div/form/textarea').send_keys(seq)
                expasy.find_element(By.XPATH, '/html/body/main/div/form/input[3]').click()
                break
            else:
                print("input box is not displayed")

    def compute():
        """get the parameters showed below"""
        while True:
            if expasy.find_element(By.XPATH, '/html/body/main').is_displayed():
                pd = {}
                pd["instability_index"] = 0
                parameters = expasy.find_element(By.XPATH, '/html/body/main').text.split("\n\n")  # 分割不同参数
                aaa = '\n'.join(parameters)
                # print(aaa)
                bbb = re.split("[:\n]", aaa)  # 将参数值 与 值分割
                noac = bbb.index("Number of amino acids") + 1
                mw = bbb.index("Molecular weight") + 1
                tp = bbb.index("Theoretical pI") + 1
                f = bbb.index("Formula") + 1
                tnoa = bbb.index("Total number of atoms") + 1
                ii = bbb.index("Instability index") + 2
                ai = bbb.index("Aliphatic index") + 1
                g = bbb.index("Grand average of hydropathicity (GRAVY)") + 1
                pd["number_of_amine_acid"] = bbb[noac].strip()
                pd["molecular_weight"] = bbb[mw].strip()
                pd["theoretical_pi"] = bbb[tp].strip()
                pd["Formula"] = bbb[f].strip()
                pd["Total_number_of_atoms"] = bbb[tnoa].strip()
                pd["instability_index"] = re.findall("\\d+\\.\\d+", bbb[ii])[0]
                pd["aliphatic_index"] = bbb[ai].strip()
                pd["gravy"] = bbb[g].strip()
                break
            else:
                print("loading")
        return pd
with open(out_file, "w", encoding='utf-8') as f:
    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        'seq_id',
        'number_of_amine_acid',
        'molecular_weight',
        'theoretical_pi',
        'Formula',
        'Total_number_of_atoms',
        'instability_index',
        'aliphatic_index',
        'gravy'))
    # 可根据自己需要添加相应的字段，
    pros = SeqIO.parse(file_name, "fasta")  # 根据自己的文件名称进行更改
    i = 0
    for pro in pros:
        print("=" * 10, "seq", i + 1, "->", pro.id, "on the way", "=" * 10)
        expasy_cal.inputseq(seq=pro.seq)
        table = expasy_cal.compute()
        number_of_amine_acid = table['number_of_amine_acid']
        molecular_weight = table['molecular_weight']
        theoretical_pi = table['theoretical_pi']
        Formula = table['Formula']
        Total_number_of_atoms = table['Total_number_of_atoms']
        instability_index = table['instability_index']
        aliphatic_index = table['aliphatic_index']
        gravy = table['gravy']
        # 对前面的字段进行赋值
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            pro.id,
            number_of_amine_acid,
            molecular_weight,
            theoretical_pi,
            Formula,
            Total_number_of_atoms,
            instability_index,
            aliphatic_index,
            gravy))  # 若需要提取的参数改变，则此处也要改变，即添加或减少{}\以及参数名称
        i += 1
        expasy.back()  # 好像back也需要重新加载页面
expasy.close()
print('完成')
