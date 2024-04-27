from TDA_fmm import Select_Best_FMM as fmm
import numpy as np

class FMM_FDR:
    def __init__(self):
        pass
        
    def get_false_model(self):
        return self.fmm.best_fmm.external_model
        
    def get_fmm_model(self):
        return self.fmm.best_fmm
    
    def fit(self, target_scores, decoy_scores, verbose = False):
        self.fmm = fmm(target_scores, decoy_scores, verbose = verbose)
        return self
        
    def estimate(self, scores):
        pep = self.fmm.best_fmm.pep(scores)
        if pep is None:
            print('target pep is none')
            return [0]*len(scores), [-1]*len(scores)
        sorted_idx = sorted(range(len(scores)), key = lambda i: scores[i], reverse = True)
        sorted_pep = pep[sorted_idx]

        # smoothing FDR into Q-value
        FDR = np.cumsum(sorted_pep) / np.arange(1,len(scores)+1)
        minFDR = 1
        for i in range(len(FDR)-1, -1, -1):
            if minFDR > FDR[i]: minFDR = FDR[i]
            else: FDR[i] = minFDR
        
        # return corresponded FDR values of input scores
        input_align_fdr = sorted(zip(sorted_idx, FDR))
        return [x[1] for x in input_align_fdr], pep

class TDA_FDR:
    def __init__(self):
        pass
    
    def fit(self, target_scores, decoy_scores, verbose = False):
        scores = list(zip(target_scores, [1]*len(target_scores)))
        scores.extend( zip(decoy_scores, [0]*len(decoy_scores)) )
        
        # target = 1, decoy = 0, keeping target before decoy
        scores.sort(reverse = True)
        
        scores = np.array(scores)
        
        target_cumsum = np.cumsum(scores[:,1])
        decoy_cumsum = np.cumsum(1-scores[:,1])
        
        fdr = decoy_cumsum / target_cumsum
        
        min_fdr = 1
        for i in range(len(fdr)-1, -1, -1):
            if fdr[i] > min_fdr:
                fdr[i] = min_fdr
            else:
                min_fdr = fdr[i]
        fdr = np.reshape(fdr, (fdr.shape[0], 1))
        self.fdr_table = np.concatenate( (scores, fdr), axis = 1)
        return self
    
    def estimate(self, scores):
        sorted_idx = sorted(range(len(scores)), key = lambda i: scores[i], reverse = True)
        FDR = [self.fdr_table[:,2].max()]*len(scores)
        i = 0;
        j = 0;
        while i < self.fdr_table.shape[0] and j < len(sorted_idx):
            if self.fdr_table[i,0] > scores[sorted_idx[j]]:
                i += 1
            elif self.fdr_table[i,0] == scores[sorted_idx[j]]:
                FDR[j] = self.fdr_table[i,2]
                j += 1
            else:
                if i == 0: FDR[j] = self.fdr_table[i,2]
                else: FDR[j] = self.fdr_table[i-1,2]
                j += 1
        
        # return corresponded FDR values of input scores with same order
        input_align_fdr = sorted(zip(sorted_idx, FDR))
        return [x[1] for x in input_align_fdr], [-1]*len(scores)
        