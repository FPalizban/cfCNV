#!/usr/bin/env python3
import csv,os
from pathlib import Path
outroot=Path(os.environ["SIM1_OUT"]); harm=Path(os.environ["SIM1_HARM"]); harm.mkdir(parents=True,exist_ok=True)
man=Path(os.environ["SIM1_MAN"])/"master.tsv"
with open(man) as f: meta={r["sample_id"]:r for r in csv.DictReader(f,delimiter="\t")}
fields=["sample_id","background","target_depth","replicate","class","tool","bin_size","chrom","start","end","cnv_type","copy_number","segment_value"]
allrows=[]
for tool in ["QDNAseq","CNVkit","ichorCNA","Spectre"]:
 rows=[]
 td=outroot/tool
 for fp in td.rglob("*.calls.tsv") if td.exists() else []:
  for r in csv.DictReader(open(fp),delimiter="\t"):
   m=meta.get(r["sample_id"],{})
   b=""
   for part in fp.parts[::-1]:
    if part.isdigit():b=part;break
   x={"sample_id":r["sample_id"],"background":m.get("background",""),"target_depth":m.get("target_depth",""),
      "replicate":m.get("replicate",""),"class":m.get("class",""),"tool":tool,"bin_size":b,
      "chrom":str(r["chrom"]).replace("chr",""),"start":int(float(r["start"])),"end":int(float(r["end"])),
      "cnv_type":r["cnv_type"].upper(),"copy_number":r.get("copy_number",""),"segment_value":r.get("segment_value","")}
   if x["cnv_type"] in ("LOSS","GAIN") and x["end"]>x["start"]:rows.append(x);allrows.append(x)
 with open(harm/f"{tool}.tsv","w",newline="") as g:
  w=csv.DictWriter(g,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
 print(tool,len(rows))
with open(harm/"all_calls.tsv","w",newline="") as g:
 w=csv.DictWriter(g,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(allrows)
print("ALL",len(allrows))
