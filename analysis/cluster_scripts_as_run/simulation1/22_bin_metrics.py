#!/usr/bin/env python3
import csv,os,re
from pathlib import Path
from collections import defaultdict
G=10000;harm=Path(os.environ["SIM1_HARM"])/"all_calls.tsv";man=Path(os.environ["SIM1_MAN"])/"master.tsv"
truthdir=Path(os.environ["SIM1_TRUTH"]);fai=Path(os.environ["REF37_FAI"]);out=Path(os.environ["SIM1_MET"])/"bin10kb_metrics.tsv"
with open(man) as f:meta={r["sample_id"]:r for r in csv.DictReader(f,delimiter="\t")}
nb=sum((int(L)+G-1)//G for c,L,*_ in (l.rstrip().split("\t") for l in open(fai)) if re.fullmatch(r"([1-9]|1[0-9]|2[0-2])",c))
def truth(rep):
 if not rep:return []
 cand=[p for p in truthdir.glob("*") if p.is_file() and re.search(rf"(rep|r)[_-]?{rep}\b",p.name,re.I)]
 if not cand:cand=list(truthdir.glob(f"*{rep}*"))
 if not cand:return []
 p=sorted(cand)[0];first=open(p).readline();delim="\t" if "\t" in first else None
 if not delim:return []
 rows=list(csv.DictReader(open(p),delimiter="\t"));ans=[]
 for x in rows:
  low={k.lower():k for k in x}
  def q(parts):
   for part in parts:
    for lk,k in low.items():
     if part in lk:return k
  ck=q(["chrom","chr"]);sk=q(["start"]);ek=q(["end"]);tk=q(["type","cnv","state","event"])
  if not all([ck,sk,ek]):continue
  t=(x.get(tk,"") if tk else "").upper()
  if any(z in t for z in ["DEL","LOSS","CN1"]):t="LOSS"
  elif any(z in t for z in ["DUP","GAIN","CN3"]):t="GAIN"
  else:continue
  ans.append((str(x[ck]).replace("chr",""),int(float(x[sk])),int(float(x[ek])),t))
 return ans
def B(iv):
 s=set()
 for c,a,b,t in iv:
  if not re.fullmatch(r"([1-9]|1[0-9]|2[0-2])",c):continue
  for k in range(a//G,(b-1)//G+1):s.add((c,k,t))
 return s
calls=defaultdict(list)
for r in csv.DictReader(open(harm),delimiter="\t"):calls[(r["sample_id"],r["tool"],r["bin_size"])].append((r["chrom"],int(r["start"]),int(r["end"]),r["cnv_type"]))
rows=[]
for (sid,tool,b),pred in calls.items():
 m=meta[sid];tb=set() if m["class"]=="negative" else B(truth(m["replicate"]));pb=B(pred)
 tp=len(tb&pb);fp=len(pb-tb);fn=len(tb-pb);ta={(c,k) for c,k,t in tb};pa={(c,k) for c,k,t in pb}
 tpa=len(ta&pa);fpa=len(pa-ta);fna=len(ta-pa);tn=max(0,nb-tpa-fpa-fna)
 pr=tp/(tp+fp) if tp+fp else (1 if not tb else 0);rec=tp/(tp+fn) if tp+fn else 1;f1=2*pr*rec/(pr+rec) if pr+rec else 0
 acc=(tpa+tn)/nb;spec=tn/(tn+fpa) if tn+fpa else 1;bal=(rec+spec)/2;den=((tpa+fpa)*(tpa+fna)*(tn+fpa)*(tn+fna))**.5;mcc=(tpa*tn-fpa*fna)/den if den else 0
 rows.append({"sample_id":sid,"background":m["background"],"target_depth":m["target_depth"],"replicate":m["replicate"],"class":m["class"],"tool":tool,"bin_size":b,
              "TP_bins":tp,"FP_bins":fp,"FN_bins":fn,"TN_bins":tn,"precision":pr,"recall":rec,"F1":f1,"accuracy":acc,"specificity":spec,"balanced_accuracy":bal,"MCC":mcc})
fields=list(rows[0]) if rows else []
with open(out,"w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
print("bin rows",len(rows))
