


#python /mnt/h/study/deep_learning/gene/project/GNNLink/OsData/MSU2RAP/MSU2RAP.py >/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/TF_os1.csv

#---------------------MSU to RAP-----------------------Convert MSU LOC-ID to RAP-DB ID: LOC_Os03g61590 -> Os03g0831400
relation={}
for i in open("/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/MSU2RAP/RAP-MSU_2018-03-29.txt"):
    rap=str(i.split()[0])
    msu=str(i.split()[1])
    if msu!="None":
        if "," in msu:
            for a in msu.split(","):
                relation[a[0:-2]] = rap
        else:
            relation[msu[0:-2]] = rap

for j in open("/mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/RNASeq/two.txt"):
    id=j.strip()
    if id in relation.keys():
        print(id,relation[id],sep="\t")
    else:
        print(id,"None",sep="\t")



#-------------------RAP to MSU-------------------------Convert RAP-DB ID to MSU LOC-ID: Os03g0831400 -> LOC_Os03g61590
# relation={}
# for i in open("/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/MSU2RAP/RAP-MSU_2018-03-29.txt"):
#         rap=str(i.split()[0])
#         msu=str(i.split()[1])
#         if rap!="None":
#             relation[rap]=msu

# for j in open("/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/deggeneid.txt"):
#     id=j.strip()
#     if id in relation.keys():
#         if "," in relation[id]:
#             s=relation[id].split(",")
#             for a in s:
#                 print(id,a,sep="\t")
#         else:
#             print(id,relation[id],sep="\t")
#     else:
#         print(id,"None",sep="\t")




