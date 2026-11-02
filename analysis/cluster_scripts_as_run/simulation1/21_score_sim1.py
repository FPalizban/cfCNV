#!/usr/bin/env python3
import csv,os,re,math
from pathlib import Path
from collections import defaultdict
harm=Path(os.environ["SIM1_HARM"])/"all_calls.tsv";truthdir=Path(os.environ["SIM1_TRUTH"])
man=Path(os.environ["SIM1_MAN"])/"master.tsv";outdir=Path(os.environ["SIM1_MET"]);outdir.mkdir(parents=True,exist_ok=True)
with open(man) as f:meta={r["sample_id"]:r for r in csv.DictReader(f,delimiter="\t")}
def truth(rep):
 if not rep:return []
 cand=[p for p in truthdir.glob("*") if p.is_file() and re.search(rf"(rep|r)[_-]?{rep}\b",p.name,re.I)]
 if not cand:cand=list(truthdir.glob(f"*{rep}*"))
 if not cand:return []
 p=sorted(cand)[0]
 with open(p) as f:first=f.readline().rstrip()
 delim="\t" if "\t" in first else None
 with open(p) as f:
  if delim:rows=list(csv.DictReader(f,delimiter="\t"))
  else:
   hdr=re.split(r"\s+",f.readline().strip());rows=[]
   for l in f:
    v=re.split(r"\s+",l.strip())
    if len(v)>=len(hdr):rows.append(dict(zip(hdr,v)))
 out=[]
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
  out.append({"chrom":str(x[ck]).replace("chr",""),"start":int(float(x[sk])),"end":int(float(x[ek])),"cnv_type":t})
 return out
def overlap(a,b):
 if a["chrom"]!=b["chrom"] or a["cnv_type"]!=b["cnv_type"]:return 0
 return max(0,min(a["end"],b["end"])-max(a["start"],b["start"]))
def recov(a,b):
 o=overlap(a,b)
 return 0 if not o else min(o/(a["end"]-a["start"]),o/(b["end"]-b["start"]))
def match(t,p,thr=.5):
 pairs=sorted([(recov(a,b),i,j) for i,a in enumerate(t) for j,b in enumerate(p) if recov(a,b)>=thr],reverse=True)
 ti=set();pj=set();m=[]
 for s,i,j in pairs:
  if i not in ti and j not in pj:ti.add(i);pj.add(j);m.append((i,j))
 return m
calls=defaultdict(list)
with open(harm) as f:
 for r in csv.DictReader(f,delimiter="\t"):
  calls[(r["sample_id"],r["tool"],r["bin_size"])].append({"chrom":r["chrom"],"start":int(r["start"]),"end":int(r["end"]),"cnv_type":r["cnv_type"]})
rows=[]
for (sid,tool,b),pred in calls.items():
 m=meta[sid];tru=[] if m["class"]=="negative" else truth(m["replicate"]);mt=match(tru,pred)
 tp=len(mt);fp=len(pred)-tp;fn=len(tru)-tp
 pr=tp/(tp+fp) if tp+fp else (1 if not tru else 0);reca=tp/(tp+fn) if tp+fn else 1;f1=2*pr*reca/(pr+reca) if pr+reca else 0
 jac=[];bp=[]
 for i,j in mt:
  t,p=tru[i],pred[j];o=overlap(t,p);u=(t["end"]-t["start"])+(p["end"]-p["start"])-o
  jac.append(o/u if u else 0);bp.append((abs(t["start"]-p["start"])+abs(t["end"]-p["end"]))/2)
 rows.append({"sample_id":sid,"background":m["background"],"target_depth":m["target_depth"],"replicate":m["replicate"],"class":m["class"],
  "tool":tool,"bin_size":b,"truth_events":len(tru),"predicted_events":len(pred),"TP":tp,"FP":fp,"FN":fn,
  "precision":pr,"recall":reca,"F1":f1,"mean_jaccard":sum(jac)/len(jac) if jac else "",
  "mean_breakpoint_error_bp":sum(bp)/len(bp) if bp else "",
  "false_altered_bp":sum(x["end"]-x["start"] for x in pred) if m["class"]=="negative" else ""})
fields=list(rows[0]) if rows else []
with open(outdir/"event_metrics.tsv","w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
print("event rows",len(rows))
