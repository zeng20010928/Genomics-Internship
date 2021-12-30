# Genomics-Internship
## 1 根据diamond结果对蛋白序列进行聚类 -- 构建基因家族
任务：构建Brevibacillus基因家族 
### 1.1 数据准备
根据链接https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/Brevibacillus
中的信息完成https://docs.qq.com/sheet/DUEZiWFBEcktGTWROzz
中表的信息
根据老师再github所给的curl命令下载本次实验所需要的数据类似下图
![image](https://user-images.githubusercontent.com/71910521/147722973-6f81c9fe-8e60-4719-a4c0-275f906e777b.png)
在下载过程中可能有些文件会出现文件下载后不能解压的情况，可能是由于文件下载不完全导致的，可以使用curl -0 --data - binary的命令来下载。 

### 1.2 对蛋白序列进行两两比对
在做diamond之前需要改一下序列名，在各自序列名后面加上GCF编号，如将`WP_003333770.1`改成`WP_003333770.1:GCF_000010165`，将所有蛋白序列合并到一个文件`all_pro.faa`。改完序列名之后的情况类似于下图
![image](https://user-images.githubusercontent.com/71910521/147723487-544bac71-c78d-4408-88d4-525591577f62.png)
改序列名有利于之后分辨获得正确的单拷贝基因家族。
* 首先先对蛋白质序列进行建库，`diamond makedb --in all_pro.faa -d allpep`，然后用diamond进行两两比对  

work_diamod.sh
```sh
#!/bin/bash
#$ -S /bin/bash
#$ -N diamond
#$ -j y
#$ -cwd
source /opt/miniconda3/bin/activate
conda activate genomelab
diamond blastp -d allpep -q /data/stdata/genomic/bioinfo2019/201941660807/shixi/data/data_80/all_pro.faa -o allBlast.tsv -f 6 #不同人数据存放位置不同，根据实际情况改路径
```
将比对的结果放入了allBlast.tsv,结果图及过程图如下
结果图
![image](https://user-images.githubusercontent.com/71910521/147726189-4816212d-ca34-4a59-966a-28d5504743a9.png)
过程图
![image](https://user-images.githubusercontent.com/71910521/147726284-aae1c773-3efc-4d32-b4d6-509603a6c2cf.png)
### 1.3 提取每个hit的score值，构建一个表征两条序列的相似性的特征值
直接使用老师所给命令
```sh
cut -f 1,2,12 allBlast.tsv > allBlast.abc
```
结果图
![image](https://user-images.githubusercontent.com/71910521/147726397-ffd63ff8-8f49-4291-9726-8bf5eb694034.png)


### 1.4 用mcl进行聚类
在集群上需要在老师的代码里加点东西
work_mcl.sh
```sh
#!/bin/bash
#$ -S /bin/bash
#$ -N MCL
#$ -j y
#$ -cwd
source /opt/miniconda3/bin/activate
conda activate genomelab
mcxload -abc allBlast.abc --stream-mirror -write-tab data.tab -o data.mci
mcl data.mci -I 1.4
mcl data.mci -I 2
mcl data.mci -I 4
mcxdump -icl out.data.mci.I14 -tabr data.tab -o dump.data.mci.I14
mcxdump -icl out.data.mci.I20 -tabr data.tab -o dump.data.mci.I20
mcxdump -icl out.data.mci.I40 -tabr data.tab -o dump.data.mci.I40
clm dist --chain out.data.mci.I{14,20,40}
```
出现三个dump文件，随机选取一个进行分析就可以了，我们组选取了dump.data.mci.I20和dump.data.mci.I40进行后续的分析。在dump文件中每一行相当于一个基因家族，我们需要筛选的单拷贝基因家族的特征为每个家族有80个成员且GCF号不能重复
据上述要求得到一个python代码如下
```python
d = "dump.data.mci.I14"#可以换成40,20
d1 = "dump.data.mci.I-14"#可以换成40,20
fo = open(d,"r")
f = open(d1,"a+")
fo = fo.readlines()
n = 0
for line in fo:
    line = line.split("\t")
    ID = []
    if len(line) == 80:
        for i in line:
            ID.append(i[15:29])
        ID_set = set(ID)
        if len(ID_set) == len(ID):
            f.write("\t".join(line))
            n += 1
f.close()
print("同学，一共{}个单拷贝基因家族".format(n))
```
根据上述代码知道dump.data.mci.I20和dump.datt.mci.I40中的单拷贝基因家族分别为72个，21个，且会生成两个包含单拷贝基因家族的文件dump.data.mci.I-20和dump.data.mci.I-40
![image](https://user-images.githubusercontent.com/71910521/147727331-a035bd02-2c0d-4941-adff-d3d9546a3712.png)


dump.data.mci.I-20，如上图一行为一个单拷贝基因家族一共72行
![image](https://user-images.githubusercontent.com/71910521/147727487-4e02e34c-1d50-4b4a-acbd-fc953cdbc267.png)


dump.data.mci.I-40，如上图
## 2.构建物种进化树
### 2.1 提取单拷贝基因家族的基因序列 
根据需求得出如下代码
```python
def extract(filename1,filename2,prefix,suffix):#提取单序列家族序列
    f = open(filename1, "r")
    f_list = f.readlines()
    f.close()
    fo = open(filename2, "r")
    fo_list = fo.readlines()
    fo.close()
    for i in range(len(f_list)):
        t = open(prefix + str(i) + suffix, "a+")
        f_name = f_list[i].strip("\n")
        f_name = f_name.split("\t")
        lines = []
        for ID in f_name:
            line = ""
            n = 0
            for row in fo_list:
                if ID in row:
                    n += 1
                    line += row
                elif row[0] != ">" and n == 1:
                    line += row
                elif row[0] == ">" and n == 1:
                    break
                else:
                    continue
            lines.append(line)
        t.write("".join(lines))
        t.close()
    return "提取完成"
#示例：获取dump.data.mci.I-40中的每个基因家族的序列
d = "dump.data.mci.I-14"#此变量可以赋值为dump.data.mci.I-20(或者40)
a = "E:/R学习、/基因组学/实习/genomelab-1/all_pro.faa"
prefix = "I14-"#可以改为I20-或者I40
suffix = ".faa"
extract(d,a,prefix,suffix)
```
会得到等于单拷贝基因家族数的序列文件。
例：I40的家族成员序列文件
![image](https://user-images.githubusercontent.com/71910521/147728152-a52ab908-4335-46d0-b86f-9d270c80189f.png)


注：为什么要选择单拷贝基因？
一般是用单拷贝的直系同源基因研究系统发生。我个人理解，不一定完全准确。
1.相对于全基因组减少计算量
2.减少旁系同源基因对序列比对和后来进化树的影响。

### 2.2 多序列比对
关于muscle的使用

使用shell脚本建立一个循环,对每个单拷贝基因家族进行多序列比对
```sh
#!/bin/bash
#$ -S /bin/bash
#$ -N muscle
#$ -j y
#$ -cwd
compack="zip"
oldsuffix="faa"
newsuffix="afa"
pack="rar"
str="muscle finshed"
dir=$(eval pwd)
unzip *.${compack}
export PATH=/data/stdata/genomic/bioinfo2019/201941632216/inter/muscle:$PATH #加入临时环境变量
for file in $(ls $dir | grep .${oldsuffix}) #构建循环
    do
        name=$(ls ${file} | cut -d. -f1) #切割获取文件名
        ./muscle -in ${file} -out ${name}.${newsuffix} #muscle运行多序列比对
	echo ${name}-${str}
	rm ${name}.${oldsuffix}
	if [ ! -f "/data/stdata/genomic/bioinfo2019/201941632216/inter/excel/${name:0:3}.${pack}" ];then #判断文件是否存在
                tar -cf ${name:0:3}.${pack} ls*.${newsuffix} #若文件不存在则创建压缩包
        else
                echo "文件存在"
        fi
	tar -crf ${name:0:3}.${pack} ls*.{newsuffix} #打包文件,将单个文件压缩进压缩包
	rm -rf ${name}.${newsuffix} #删除多序列比对多余结果
    done
echo ${str}
```
I20的多序列比对后的结果文件数量
![image](https://user-images.githubusercontent.com/71910521/147729055-6cad902c-5e80-4c52-bf6e-3c686bf41606.png)


文件内容（只取了其中一个），经过多序列比对之后会将不等长的序列使用gap填充，使之等长。
![image](https://user-images.githubusercontent.com/71910521/147729214-df61789e-fd85-45b7-bd8f-2c0c13d35fee.png)


I40的多序列比对后的结果文件数量
![image](https://user-images.githubusercontent.com/71910521/147730200-54c0769b-6a17-467a-ad6b-10733fb5bf98.png)


文件内容
![image](https://user-images.githubusercontent.com/71910521/147730315-f1fef09e-0a2f-402e-a25c-4f76ef782fbf.png)
### 2.3 比对后序列合并  
对于每个菌株的ID（GCF号），都需要在每一个muscle比对文件中搜索对应的序列，将序列逐个粘贴成1条完整的序列
如果我们有80个菌株，那么合并的文件应该有160行（fasta格式，一行ID、一行序列）
由上述要求所得代码如下
```python
def filenames(path,prefix,suffix):
    import os
    path_list = os.listdir(path)
    path_list1 = []
    for f in path_list:
        if f.find(suffix) != -1:
            path_list1.append(f)
    filename = []
    for i in range(len(path_list1)):
        filename.append(prefix + str(i) + suffix)
    return filename
def ID(filename,IDname):
    f = open(filename,"r").readlines()
    fo = open(IDname,"a+")
    for line in f:
        if line[0] == ">":
            fo.write(">"+line[16:29]+"\n")
    fo.close()
    return
def merge_seq(ID_merge_file,new_file,total_mergefilename):
    f = open(new_file, "a+")
    fo = open(ID_merge_file, "r")
    fo_list = fo.readlines()
    for i in fo_list:
        lines = i
        for I in total_mergefilename:
            seq = open(I, "r").readlines()
            n = 0
            for row in seq:
                if row.find(i[1:16].strip("\n")) != -1:
                    n += 1
                elif row[0] != ">" and n == 1:
                    row = row.strip("\n")
                    lines += row
                elif row[0] == ">" and n == 1:
                    break
                else:
                    continue
        f.write(lines + "\n")
    f.close()
    return "序列融合完毕"
if __name__ == '__main__':
    path = "E:/R学习、/基因组学/实习/results/I14"  # 20,40
    prefix = "I14-"  # 20,40
    suffix = ".afa"
    file = "I14-0.afa"  # 从中提取80个GCP名
    name = "ID_I14.faa"  # 80个GCP名字所在文件名字，此时为空
    I_filename = filenames(path, prefix, suffix)  # 获得需要合并的文件的总名称数且进行了顺序排列
    ID(file, name)  # 将80个GCP的ID储存在了name代表的文件名中
    after_mergefile = "all_I14.faa"  # 序列融合后得到的文件名称
    merge_seq(name, after_mergefile, I_filename)
```
分别运行完上述代码后得到两种聚类分别合并后的文件all_I20.faa和all_I40.faa
文件内容如图
all_I20.faa（部分）
![image](https://user-images.githubusercontent.com/71910521/147731455-833040a0-bdff-4a4f-a6b4-b487a250d21e.png)
all_I40.faa（部分）
![image](https://user-images.githubusercontent.com/71910521/147731597-c7ffb052-8f63-4909-8172-b3ac5abd98d9.png)

### 2.4 构建进化树
集群上有FastTree软件，可直接使用fastree protein_alignment > tree
最后，将tree文件放到mega或者其他的构树软件中，出图如下
I20
![0Z~KO)A0O2RAIHRTTM4(FVA](https://user-images.githubusercontent.com/71910521/147730666-5eb0ee56-c738-4ac0-9372-93e33a8e51df.png)

I40
![TL}ZGZ9KPS9 FT06BWNQV(E](https://user-images.githubusercontent.com/71910521/147730695-d570001b-83f8-4c15-96c5-6bc924dd1ecb.png)
