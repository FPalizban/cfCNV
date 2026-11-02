#!/usr/bin/env python3
import argparse,csv,math,os,re,statistics,sys
from pathlib import Path
from collections import defaultdict


def norm(s): return re.sub(r'[^a-z0-9]+','',str(s).lower())
def sf(x):
    try:
        if x is None or str(x).strip() in ('','.','NA','NaN','nan'): return None
        return float(str(x).strip())
    except: return None
def si(x):
    v=sf(x); return int(round(v)) if v is not None else None
def div(a,b): return None if b in (None,0) else a/b
def f1(p,r):
    if p is None or r is None: return None
    return 0.0 if p+r==0 else 2*p*r/(p+r)
def mcc(tp,tn,fp,fn):
    d=(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
    return None if d<=0 else (tp*tn-fp*fn)/math.sqrt(d)
def med(xs):
    xs=[x for x in xs if x is not None]; return statistics.median(xs) if xs else None
def cchr(x):
    s=re.sub(r'^chr','',str(x).strip(),flags=re.I)
    return s if s in {str(i) for i in range(1,23)} else None
def ctype(v):
    if v is None: return None
    s=str(v).strip().lower()
    if any(k in s for k in ('del','loss','deletion')): return 'loss'
    if any(k in s for k in ('dup','gain','amp','amplification')): return 'gain'
    v=sf(s)
    if v is not None: return 'loss' if v<2 else ('gain' if v>2 else 'neutral')
    return None

def read_table(path):
    with open(path,newline='') as f:
        first=f.readline()
        if not first: return [],[]
        d='\t' if first.count('\t')>=first.count(',') else ','; f.seek(0)
        r=csv.DictReader(f,delimiter=d); return r.fieldnames or [],list(r)
def pick(fields,cands):
    m={norm(x):x for x in fields}
    for c in cands:
        if norm(c) in m: return m[norm(c)]
    for c in cands:
        nc=norm(c)
        for k,v in m.items():
            if nc and (k.endswith(nc) or nc in k): return v
    return None

def sample_meta(path):
    f,rows=read_table(path)
    sc=pick(f,['sample_id','sample','id']); cc=pick(f,['class','sample_class','status','label'])
    dc=pick(f,['depth_x','depth','target_depth']); bc=pick(f,['background','source_sample','source'])
    rc=pick(f,['replicate','rep','placement_replicate']); tc=pick(f,['truth_path','truth_bed','truth','ground_truth'])
    if not sc: raise RuntimeError(f'No sample column in {path}: {f}')
    out={}
    for r in rows:
        s=str(r.get(sc,'')).strip()
        if not s: continue
        raw=str(r.get(cc,'')).lower() if cc else ''
        cl='negative' if ('neg' in raw or 'control' in raw or (not raw and re.search('neg|control',s,re.I))) else 'positive'
        d=sf(r.get(dc)) if dc else None
        if d is None:
            m=re.search(r'(\d+(?:\.\d+)?)x(?:\b|_)',s,re.I); d=float(m.group(1)) if m else None
        bg=str(r.get(bc,'')).strip() if bc else ''
        if not bg:
            m=re.search(r'(IH0[123])',s,re.I); bg=m.group(1).upper() if m else ''
        rep=si(r.get(rc)) if rc else None
        if rep is None:
            m=re.search(r'(?:rep(?:licate)?|r)[_-]?([1-5])',s,re.I); rep=int(m.group(1)) if m else None
        out[s]={'class':cl,'depth_x':d,'background':bg,'replicate':rep,'truth_path':str(r.get(tc,'')).strip() if tc else ''}
    return out

def pbin(v):
    if v is None: return None
    s=str(v).strip().lower().replace(',','')
    m=re.fullmatch(r'(\d+(?:\.\d+)?)\s*(kb|k|mb|m|bp)?',s)
    if not m: return si(s)
    x=float(m.group(1)); u=m.group(2) or ''
    if u in ('kb','k'): x*=1e3
    elif u in ('mb','m'): x*=1e6
    return int(round(x))

def load_tasks(tool,path):
    f,rows=read_table(path); sc=pick(f,['sample_id','sample','id']) or f[0]
    bc=pick(f,['bin_size_bp','bin_bp','bin','binsize','target_avg_size','window','window_size'])
    out=[]
    for i,r in enumerate(rows):
        s=str(r.get(sc,'')).strip(); b=pbin(r.get(bc)) if bc else None
        if tool in ('QDNAseq','CNVkit') and b is not None and b<=2000: b*=1000
        if tool=='ichorCNA' and b is None: b=1_000_000
        if tool=='Spectre' and b is None: b=1_000
        if s: out.append({'tool':tool,'sample':s,'bin_bp':b,'task_index':i})
    return out

def parse_bedlike(path,explicit=False):
    with open(path,errors='replace') as f: lines=[x.rstrip('\n') for x in f if x.strip() and not x.startswith('##')]
    if not lines: return []
    if any(x.startswith('#CHROM') for x in lines) or str(path).lower().endswith('.vcf'):
        out=[]
        for ln in lines:
            if ln.startswith('#'): continue
            p=ln.split('\t')
            if len(p)<8: continue
            ch=cchr(p[0]); pos=si(p[1]); info={}
            for kv in p[7].split(';'):
                if '=' in kv:
                    k,v=kv.split('=',1); info[k.upper()]=v
            en=si(info.get('END')); ty=ctype(info.get('SVTYPE')) or ctype(info.get('CN'))
            if ch and pos and en and en>pos and ty in ('gain','loss'):
                out.append({'chrom':ch,'start':pos-1,'end':en,'type':ty})
        return out
    first=lines[0]; d='\t' if first.count('\t')>=first.count(',') else ','
    toks=[norm(x) for x in first.split(d)]
    header=any(x in toks for x in ('chrom','chromosome','chr','start','end','stop','locstart','locend','svtype','type','call','cn','copynumber','segmean','log2'))
    out=[]
    if header:
        r=csv.DictReader(lines,delimiter=d); flds=r.fieldnames or []
        ch=pick(flds,['chromosome','chrom','chr']); st=pick(flds,['start','loc.start','begin']); en=pick(flds,['end','stop','loc.end'])
        ty=pick(flds,['svtype','type','call','state','event_type','cnv_type']); cn=pick(flds,['copy.number','copy_number','copynumber','cn','corrected_copy_number'])
        lg=pick(flds,['seg.mean','segmean','log2','log2ratio','segment_mean','logr','mean'])
        n_valid_intervals=0
        for z in r:
            C=cchr(z.get(ch,'')) if ch else None; S=si(z.get(st)) if st else None; E=si(z.get(en)) if en else None
            if not C or S is None or E is None or E<=S: continue
            n_valid_intervals += 1
            T=ctype(z.get(ty)) if ty else None
            if T is None and cn: T=ctype(z.get(cn))
            if T is None and lg and explicit:
                v=sf(z.get(lg)); T='gain' if v is not None and v>0 else ('loss' if v is not None and v<0 else None)
            if T in ('gain','loss'): out.append({'chrom':C,'start':S,'end':E,'type':T})
        if n_valid_intervals>0 and not explicit and ty is None and cn is None:
            raise RuntimeError('Segment table has intervals but no explicit CNV type/copy-number column; refusing to interpret raw segments as calls')
        return out
    for ln in lines:
        p=ln.split(d)
        if len(p)<3: continue
        C=cchr(p[0]); S=si(p[1]); E=si(p[2]); T=None
        for x in p[3:]:
            sx=str(x).strip().lower()
            if any(k in sx for k in ('del','loss','deletion')): T='loss'; break
            if any(k in sx for k in ('dup','gain','amp','amplification')): T='gain'; break
        if C and S is not None and E is not None and E>S and T: out.append({'chrom':C,'start':S,'end':E,'type':T})
    return out

def truth_index(meta,truth_dir):
    files=[str(p) for p in Path(truth_dir).rglob('*') if p.is_file()]
    byrep=defaultdict(list)
    for p in files:
        m=re.search(r'(?:rep(?:licate)?|r)[_-]?([1-5])',os.path.basename(p),re.I) or re.search(r'truth[_-]?([1-5])',os.path.basename(p),re.I)
        if m: byrep[int(m.group(1))].append(p)
    cache={}; out={}; src={}
    def get(p):
        if p not in cache: cache[p]=parse_bedlike(p,True)
        return cache[p]
    for s,m in meta.items():
        if m['class']=='negative': out[s]=[]; src[s]=''; continue
        p=m['truth_path']
        if p and os.path.exists(p): out[s]=get(p); src[s]=p; continue
        cand=sorted(byrep.get(m['replicate'],[]),key=lambda x:(0 if 'truth' in os.path.basename(x).lower() else 1,len(x)))
        good=None
        for c in cand:
            try:
                if get(c): good=c; break
            except: pass
        out[s]=get(good) if good else None; src[s]=good or ''
    return out,src

ALIASES={'QDNAseq':['QDNAseq','qdnaseq'],'CNVkit':['CNVkit','cnvkit'],'ichorCNA':['ichorCNA','ichorcna','ichor'],'Spectre':['Spectre','spectre']}
def root_for(base,tool):
    for a in ALIASES[tool]:
        p=Path(base)/a
        if p.exists(): return p
    return Path(base)/ALIASES[tool][0]
def files_for(root,s,b):
    ds=[]
    if b is not None:
        for q in (str(b),str(int(b/1000)),f'{int(b/1000)}kb',f'{b}bp'):
            p=root/s/q
            if p.exists(): ds.append(p)
    p=root/s
    if p.exists(): ds.append(p)
    out=[]; seen=set()
    for d in ds:
        for f in d.rglob('*'):
            if f.is_file() and f not in seen: seen.add(f); out.append(f)
    return out
def choose(tool,fs):
    def rank(p):
        n=p.name.lower()
        if any(n.endswith(x) for x in ('.png','.pdf','.jpg','.rds','.rdata','.log','.err','.out')): return 99
        if tool=='QDNAseq':
            if n.endswith('.calls.vcf'): return 0
            if n.endswith('.calls.tsv'): return 1
        if tool=='CNVkit':
            if 'call' in n and n.endswith('.cns'): return 0
            if n.endswith('.calls.tsv'): return 1
            if n.endswith('.cns'): return 4
        if tool=='ichorCNA':
            if 'seg' in n and n.endswith(('.txt','.tsv','.seg')): return 0
            if 'cna' in n and n.endswith(('.txt','.tsv')): return 1
        if tool=='Spectre':
            if n.endswith('.vcf') and 'cnv' in n: return 0
            if n.endswith('.vcf'): return 1
        return 10 if n.endswith(('.vcf','.bed','.tsv','.txt','.csv')) else 50
    fs=sorted(fs,key=lambda p:(rank(p),len(str(p))))
    return fs[0] if fs and rank(fs[0])<50 else None

def ov(a,b): return 0 if a['chrom']!=b['chrom'] else max(0,min(a['end'],b['end'])-max(a['start'],b['start']))
def match(T,C,thr=.5,sametype=True):
    cand=[]
    for i,t in enumerate(T):
        lt=t['end']-t['start']
        for j,c in enumerate(C):
            if t['chrom']!=c['chrom'] or (sametype and t['type']!=c['type']): continue
            o=ov(t,c)
            if not o: continue
            lc=c['end']-c['start']; rt=o/lt; rc=o/lc
            if rt>=thr and rc>=thr: cand.append((o,min(rt,rc),i,j,rt,rc))
    cand.sort(reverse=True); ti=set(); ci=set(); out=[]
    for o,_,i,j,rt,rc in cand:
        if i in ti or j in ci: continue
        ti.add(i); ci.add(j); t=T[i]; c=C[j]
        out.append({'ti':i,'ci':j,'ov':o,'rt':rt,'rc':rc,'se':abs(c['start']-t['start']),'ee':abs(c['end']-t['end']),
                    'be':(abs(c['start']-t['start'])+abs(c['end']-t['end']))/2,'ts':t['end']-t['start'],'cs':c['end']-c['start'],
                    'sze':abs((c['end']-c['start'])-(t['end']-t['start'])),'rsze':abs((c['end']-c['start'])-(t['end']-t['start']))/(t['end']-t['start']),
                    'typeok':t['type']==c['type']})
    return out

def merge(X):
    d=defaultdict(list)
    for x in X: d[x['chrom']].append((x['start'],x['end']))
    out={}
    for ch,xs in d.items():
        xs.sort(); m=[]
        for s,e in xs:
            if not m or s>m[-1][1]: m.append([s,e])
            else: m[-1][1]=max(m[-1][1],e)
        out[ch]=m
    return out
def tlen(m): return sum(e-s for xs in m.values() for s,e in xs)
def ilen(a,b):
    z=0
    for ch in set(a)&set(b):
        x,y=a[ch],b[ch]; i=j=0
        while i<len(x) and j<len(y):
            s=max(x[i][0],y[j][0]); e=min(x[i][1],y[j][1]); z+=max(0,e-s)
            if x[i][1]<y[j][1]: i+=1
            else: j+=1
    return z
def bpm(T,C,G):
    a,b=merge(T),merge(C); tb,cb=tlen(a),tlen(b); tp=ilen(a,b); fn=tb-tp; fp=cb-tp; union=tb+cb-tp; tn=max(0,G-union)
    p=div(tp,tp+fp); r=div(tp,tp+fn); sp=div(tn,tn+fp)
    return {'truth_bp':tb,'called_bp':cb,'overlap_bp':tp,'bp_TP':tp,'bp_FP':fp,'bp_TN':tn,'bp_FN':fn,'bp_precision':p,'bp_recall':r,
            'bp_f1':f1(p,r),'bp_specificity':sp,'bp_accuracy':div(tp+tn,G),'bp_balanced_accuracy':(r+sp)/2 if r is not None and sp is not None else None,
            'bp_MCC':mcc(tp,tn,fp,fn),'bp_jaccard':div(tp,union)}
def grid_ranges(X,lengths,B=10000,frac=.5):
    """Return per-chrom inclusive ranges of bins covered by >=frac of that bin."""
    M=merge(X); ranges=defaultdict(list); partial=defaultdict(int)
    for ch,xs in M.items():
        L=lengths.get(ch)
        if not L: continue
        for s,e in xs:
            s=max(0,s); e=min(L,e)
            if e<=s: continue
            b0=s//B; b1=(e-1)//B
            if b0==b1:
                partial[(ch,b0)] += e-s
                continue
            cov0=(b0+1)*B-s; len0=min(B,L-b0*B)
            if cov0 >= frac*len0: ranges[ch].append([b0,b0])
            cov1=e-b1*B; len1=min(B,L-b1*B)
            if cov1 >= frac*len1: ranges[ch].append([b1,b1])
            if b1-b0>1: ranges[ch].append([b0+1,b1-1])
    for (ch,b),cov in partial.items():
        L=lengths[ch]; bl=min(B,L-b*B)
        if cov >= frac*bl: ranges[ch].append([b,b])
    out={}
    for ch,rs in ranges.items():
        rs.sort(); mm=[]
        for a,b in rs:
            if not mm or a>mm[-1][1]+1: mm.append([a,b])
            else: mm[-1][1]=max(mm[-1][1],b)
        out[ch]=mm
    return out

def range_count(R): return sum(b-a+1 for rs in R.values() for a,b in rs)
def range_intersection(A,B):
    n=0
    for ch in set(A)&set(B):
        a,b=A[ch],B[ch]; i=j=0
        while i<len(a) and j<len(b):
            lo=max(a[i][0],b[j][0]); hi=min(a[i][1],b[j][1])
            if hi>=lo: n+=hi-lo+1
            if a[i][1]<b[j][1]: i+=1
            else: j+=1
    return n

def grid_metrics(T,C,lengths,B=10000,frac=.5):
    A=grid_ranges(T,lengths,B,frac); D=grid_ranges(C,lengths,B,frac)
    tb=range_count(A); cb=range_count(D); tp=range_intersection(A,D); fn=tb-tp; fp=cb-tp
    total=sum((L+B-1)//B for L in lengths.values()); union=tb+cb-tp; tn=max(0,total-union)
    p=div(tp,tp+fp); r=div(tp,tp+fn); sp=div(tn,tn+fp)
    return {'grid10k_truth_positive_bins':tb,'grid10k_called_positive_bins':cb,'grid10k_TP':tp,'grid10k_FP':fp,'grid10k_TN':tn,'grid10k_FN':fn,
            'grid10k_precision':p,'grid10k_recall':r,'grid10k_f1':f1(p,r),'grid10k_specificity':sp,'grid10k_accuracy':div(tp+tn,total),
            'grid10k_balanced_accuracy':(r+sp)/2 if r is not None and sp is not None else None,'grid10k_MCC':mcc(tp,tn,fp,fn)}

def write_csv(p,rows,cols):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore'); w.writeheader()
        for r in rows: w.writerow({k:'' if r.get(k) is None else r.get(k) for k in cols})

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--sim1-out',required=True); ap.add_argument('--sim1-man',required=True); ap.add_argument('--truth-dir',required=True); ap.add_argument('--fai',required=True); ap.add_argument('--outdir',required=True); a=ap.parse_args()
    meta=sample_meta(Path(a.sim1_man)/'master.tsv')
    tasks=[]
    for tool,name in [('QDNAseq','qdnaseq_tasks.tsv'),('CNVkit','cnvkit_tasks.tsv'),('ichorCNA','ichor_tasks.tsv'),('Spectre','spectre_tasks.tsv')]:
        x=load_tasks(tool,Path(a.sim1_man)/name); tasks+=x; print(tool,'expected',len(x),file=sys.stderr)
    T,Ts=truth_index(meta,a.truth_dir)
    badtruth=[s for s,m in meta.items() if m['class']=='positive' and T.get(s) is None]
    if badtruth:
        Path(a.outdir).mkdir(parents=True,exist_ok=True); (Path(a.outdir)/'unresolved_truth_samples.txt').write_text('\n'.join(badtruth)+'\n'); raise RuntimeError('Unresolved truth samples; see unresolved_truth_samples.txt')
    bad_design=[]
    for ss,mm in meta.items():
        if mm['class']!='positive': continue
        tt=T[ss]
        ng=sum(x['type']=='gain' for x in tt); nl=sum(x['type']=='loss' for x in tt)
        if len(tt)!=30 or ng!=15 or nl!=15:
            bad_design.append(f'{ss}\tN={len(tt)}\tgain={ng}\tloss={nl}\t{Ts.get(ss,"")}')
    if bad_design:
        Path(a.outdir).mkdir(parents=True,exist_ok=True); (Path(a.outdir)/'truth_design_mismatch.txt').write_text('\n'.join(bad_design)+'\n')
        raise RuntimeError('Truth parser did not recover the frozen 30-event (15 gain/15 loss) design; see truth_design_mismatch.txt')
    gl={}
    with open(a.fai) as f:
        for ln in f:
            p=ln.split('\t'); ch=cchr(p[0]);
            if ch: gl[ch]=int(p[1])
    G=sum(gl[str(i)] for i in range(1,23)); Ggb=G/1e9
    rows=[]; status=[]; callsout=[]; matchesout=[]; detectout=[]
    for n,q in enumerate(tasks,1):
        tool,s,b=q['tool'],q['sample'],q['bin_bp']; m=meta[s]; truth=T[s] or []; root=root_for(a.sim1_out,tool); fn=choose(tool,files_for(root,s,b))
        if fn is None: calls=None; st='MISSING_OUTPUT'; msg='No parseable call file found'
        else:
            try:
                calls=parse_bedlike(fn,explicit=('call' in fn.name.lower() or fn.suffix.lower()=='.vcf' or tool=='Spectre')); st='OK_ZERO_CALLS' if len(calls)==0 else 'OK'; msg=''
            except Exception as e: calls=None; st='PARSE_ERROR'; msg=repr(e)
        status.append({'tool':tool,'sample':s,'analysis_bin_bp':b,'task_index':q['task_index'],'status':st,'output_file':str(fn) if fn else '','message':msg})
        r={'tool':tool,'sample':s,'background':m['background'],'sample_class':m['class'],'replicate':m['replicate'],'depth_x':m['depth_x'],'analysis_bin_bp':b,'analysis_bin_kb':b/1000 if b else None,
           'task_index':q['task_index'],'output_status':st,'output_file':str(fn) if fn else '','truth_file':Ts[s],'truth_events':len(truth),'truth_gains':sum(x['type']=='gain' for x in truth),'truth_losses':sum(x['type']=='loss' for x in truth)}
        if calls is None: rows.append(r); continue
        r.update({'called_events':len(calls),'called_gains':sum(x['type']=='gain' for x in calls),'called_losses':sum(x['type']=='loss' for x in calls)})
        for j,c in enumerate(calls): callsout.append({'tool':tool,'sample':s,'background':m['background'],'sample_class':m['class'],'replicate':m['replicate'],'depth_x':m['depth_x'],'analysis_bin_bp':b,'call_index':j,**c,'size_bp':c['end']-c['start'],'output_file':str(fn)})
        ms={th:match(truth,calls,th,True) for th in (.1,.3,.5)}
        for th in ms: r[f'overlap_events_{int(th*100)}pct']=len(ms[th])
        detected={th:{z['ti']:z for z in ms[th]} for th in ms}
        for ti,t in enumerate(truth):
            z50=detected[.5].get(ti)
            detectout.append({'tool':tool,'sample':s,'background':m['background'],'sample_class':m['class'],'replicate':m['replicate'],'depth_x':m['depth_x'],'analysis_bin_bp':b,
                'truth_index':ti,'chrom':t['chrom'],'truth_start':t['start'],'truth_end':t['end'],'truth_size_bp':t['end']-t['start'],'truth_type':t['type'],
                'detected_10pct':int(ti in detected[.1]),'detected_30pct':int(ti in detected[.3]),'detected_50pct':int(ti in detected[.5]),
                'matched_call_index_50pct':z50['ci'] if z50 else None,'overlap_bp_50pct':z50['ov'] if z50 else None,'breakpoint_error_bp_50pct':z50['be'] if z50 else None,'relative_size_error_50pct':z50['rsze'] if z50 else None})
        mm=ms[.5]; tp=len(mm); fp=len(calls)-tp; fnn=len(truth)-tp; p=div(tp,tp+fp); rec=div(tp,tp+fnn)
        r.update({'TP_events_50pct':tp,'FP_events_50pct':fp,'FN_events_50pct':fnn,'precision_50pct':p,'recall_50pct':rec,'F1_50pct':f1(p,rec),'FDR_50pct':div(fp,tp+fp),
                  'truth_detected_fraction_50pct':div(tp,len(truth)),'false_calls':fp,'false_calls_per_Gb':fp/Ggb,'median_breakpoint_error_bp':med([x['be'] for x in mm]),
                  'median_start_error_bp':med([x['se'] for x in mm]),'median_end_error_bp':med([x['ee'] for x in mm]),'median_size_error_bp':med([x['sze'] for x in mm]),'median_relative_size_error':med([x['rsze'] for x in mm])})
        cm=match(truth,calls,.5,False); tc=sum(x['typeok'] for x in cm); r.update({'coordinate_matches_50pct':len(cm),'type_correct_matches_50pct':tc,'gain_loss_accuracy_50pct':div(tc,len(cm))}); r.update(bpm(truth,calls,G)); r.update(grid_metrics(truth,calls,gl,10000,.5)); r['false_altered_fraction']=div(r.get('bp_FP'),G); r['called_genome_fraction']=div(r.get('called_bp'),G)
        for z in mm:
            t=truth[z['ti']]; c=calls[z['ci']]; matchesout.append({'tool':tool,'sample':s,'background':m['background'],'sample_class':m['class'],'replicate':m['replicate'],'depth_x':m['depth_x'],'analysis_bin_bp':b,'threshold':.5,
                'truth_index':z['ti'],'call_index':z['ci'],'chrom':t['chrom'],'truth_start':t['start'],'truth_end':t['end'],'truth_type':t['type'],'truth_size_bp':z['ts'],'call_start':c['start'],'call_end':c['end'],'call_type':c['type'],'call_size_bp':z['cs'],'overlap_bp':z['ov'],'truth_overlap_fraction':z['rt'],'call_overlap_fraction':z['rc'],'start_error_bp':z['se'],'end_error_bp':z['ee'],'breakpoint_error_bp':z['be'],'size_error_bp':z['sze'],'relative_size_error':z['rsze']})
        rows.append(r)
        if n%100==0: print('processed',n,'/',len(tasks),file=sys.stderr)
    cols=['tool','sample','background','sample_class','replicate','depth_x','analysis_bin_bp','analysis_bin_kb','task_index','output_status','output_file','truth_file','truth_events','truth_gains','truth_losses','called_events','called_gains','called_losses','overlap_events_10pct','overlap_events_30pct','overlap_events_50pct','TP_events_50pct','FP_events_50pct','FN_events_50pct','precision_50pct','recall_50pct','F1_50pct','FDR_50pct','truth_detected_fraction_50pct','false_calls','false_calls_per_Gb','coordinate_matches_50pct','type_correct_matches_50pct','gain_loss_accuracy_50pct','median_breakpoint_error_bp','median_start_error_bp','median_end_error_bp','median_size_error_bp','median_relative_size_error','truth_bp','called_bp','overlap_bp','bp_TP','bp_FP','bp_TN','bp_FN','bp_precision','bp_recall','bp_f1','bp_specificity','bp_accuracy','bp_balanced_accuracy','bp_MCC','bp_jaccard','grid10k_truth_positive_bins','grid10k_called_positive_bins','grid10k_TP','grid10k_FP','grid10k_TN','grid10k_FN','grid10k_precision','grid10k_recall','grid10k_f1','grid10k_specificity','grid10k_accuracy','grid10k_balanced_accuracy','grid10k_MCC','false_altered_fraction','called_genome_fraction']
    out=Path(a.outdir); write_csv(out/'simulation1_overlap_master.csv',rows,cols)
    write_csv(out/'run_status.csv',status,['tool','sample','analysis_bin_bp','task_index','status','output_file','message'])
    write_csv(out/'harmonized_calls.csv',callsout,['tool','sample','background','sample_class','replicate','depth_x','analysis_bin_bp','call_index','chrom','start','end','size_bp','type','output_file'])
    write_csv(out/'event_matches_50pct.csv',matchesout,['tool','sample','background','sample_class','replicate','depth_x','analysis_bin_bp','threshold','truth_index','call_index','chrom','truth_start','truth_end','truth_type','truth_size_bp','call_start','call_end','call_type','call_size_bp','overlap_bp','truth_overlap_fraction','call_overlap_fraction','start_error_bp','end_error_bp','breakpoint_error_bp','size_error_bp','relative_size_error'])
    write_csv(out/'truth_event_detection.csv',detectout,['tool','sample','background','sample_class','replicate','depth_x','analysis_bin_bp','truth_index','chrom','truth_start','truth_end','truth_size_bp','truth_type','detected_10pct','detected_30pct','detected_50pct','matched_call_index_50pct','overlap_bp_50pct','breakpoint_error_bp_50pct','relative_size_error_50pct'])
    neg=[x for x in rows if x.get('sample_class')=='negative']
    write_csv(out/'negative_control_metrics.csv',neg,cols)
    valid=[x for x in rows if x.get('output_status') in ('OK','OK_ZERO_CALLS')]
    summary=[]
    for tool in ('QDNAseq','CNVkit','ichorCNA','Spectre'):
        rr=[x for x in valid if x['tool']==tool]; pp=[x for x in rr if x['sample_class']=='positive']; nn=[x for x in rr if x['sample_class']=='negative']
        summary.append({'tool':tool,'valid_configurations':len(rr),'positive_configurations':len(pp),'negative_configurations':len(nn),
            'median_precision_50pct':med([x.get('precision_50pct') for x in pp]),'median_recall_50pct':med([x.get('recall_50pct') for x in pp]),'median_F1_50pct':med([x.get('F1_50pct') for x in pp]),
            'median_bp_balanced_accuracy':med([x.get('bp_balanced_accuracy') for x in rr]),'median_bp_MCC':med([x.get('bp_MCC') for x in rr]),'median_bp_jaccard':med([x.get('bp_jaccard') for x in rr]),
            'median_breakpoint_error_bp':med([x.get('median_breakpoint_error_bp') for x in pp]),'median_negative_false_calls_per_Gb':med([x.get('false_calls_per_Gb') for x in nn])})
    write_csv(out/'tool_summary.csv',summary,list(summary[0].keys()))
    cnt=defaultdict(int)
    for x in status: cnt[(x['tool'],x['status'])]+=1
    with open(out/'SCORING_SUMMARY.txt','w') as f:
        f.write(f'Master samples: {len(meta)}\nExpected configurations: {len(tasks)}\nAutosomal genome bp: {G}\n\n')
        for k in sorted(cnt): f.write(f'{k[0]}\t{k[1]}\t{cnt[k]}\n')
        f.write('\nPrimary event criterion: one-to-one, same chromosome, same gain/loss direction, >=50% reciprocal overlap.\n10% and 30% sensitivity analyses are also included.\nMissing outputs are not treated as zero calls.\n')
    bad=[x for x in status if x['status'] not in ('OK','OK_ZERO_CALLS')]
    if bad:
        print(f'WARNING: {len(bad)} missing/unparseable configurations; see run_status.csv',file=sys.stderr); sys.exit(3)
    print('DONE',out,file=sys.stderr)
if __name__=='__main__': main()
