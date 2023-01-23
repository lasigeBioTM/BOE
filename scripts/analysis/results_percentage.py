# This script takes the results from the output files and verifies the gold standard.
# Logs percentage of relations identified in abstracts.
import sys

def get_true_results(gs, results, outfile):
    gs_rl_d = {}
    gs_d = {}
    gs_lines = gs.readlines()
    gs.close()
    for line in gs_lines:
        line = line.strip()
        line = line.split('\t')
        if line[0] in gs_d:
            gs_d[line[0]].append((line[1],line[2],line[3]))
            if line[1] in gs_rl_d[line[0]]:
                gs_rl_d[line[0]][line[1]] += 1
            else:
                gs_rl_d[line[0]][line[1]] = 1
        else:
            gs_rl_d[line[0]] = {line[1] : 1}
            gs_d[line[0]] = [(line[1],line[2],line[3])]

    p_rl_d = {}
    p_d = {}
    results_lines = results.readlines()
    results.close()
    for line in results_lines:
        line = line.strip()
        line = line.split('\t')
        if line[0] in gs_d and (line[1],line[2],line[3]) in gs_d[line[0]]:
            if line[0] in p_d:
                p_d[line[0]].append((line[1],line[2],line[3]))
                if line[1] in p_rl_d[line[0]]:
                    p_rl_d[line[0]][line[1]] += 1
                else:
                    p_rl_d[line[0]][line[1]] = 1
            else:
                p_rl_d[line[0]] = {line[1] : 1}
                p_d[line[0]] = [(line[1],line[2],line[3])]
    
    lfp = 0 #lfp - lower than fifty percent
    hfp = 0 #hfp - higher or equal than fifty percent
    ohp = 0 #hp - one hundred percent
    for abs, rels in p_d.items():
        # pp - percentage points
        pp = (sum(p_rl_d[abs].values())/sum(gs_rl_d[abs].values()))*100
        if pp >= 50 and pp < 100:
            hfp += 1
        elif pp < 50:
            lfp += 1
        elif pp == 100:
            ohp += 1
        outfile.write('%s\t%d%%\n' % 
                        (abs, pp))
        
        for rel in rels:
            outfile.write('%s\t%s\t%s\t%s\n' %
                            (abs, rel[0], rel[1], rel[2]))
    print(lfp, hfp, ohp)

gs = open(sys.argv[1], 'r', encoding='utf-8')
results = open(sys.argv[2], 'r', encoding='utf-8')
outfile = open(sys.argv[3], 'w', encoding='utf-8')
get_true_results(gs, results, outfile)