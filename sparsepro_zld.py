import pandas as pd
import argparse
import os
import numpy as np
from scipy.special import softmax
from scipy.stats import chi2
from scipy.stats import chi2_contingency
from scipy.stats import entropy

np.set_printoptions(precision=4, linewidth=200)


def title():
    print('**********************************************************************')
    print('* SparsePro for efficient locus-wide fine-mapping                    *')
    print('* Version 3.0.0                                                      *')
    print('* (C) Wenmin Zhang (wenmin.zhang@mail.mcgill.ca)                     *')
    print('**********************************************************************')
    print()


class SparsePro(object):

    def __init__(self, P, K, XX, h2, var_b):
        """initialize and set hyperparameters"""
        self.p = P
        self.k = K
        self.gamma = np.zeros((self.p, self.k))
        self.beta_mu = np.zeros((self.p, self.k))
        self.beta_prior_tau = np.tile((1.0 / var_b * np.array([k + 1 for k in range(self.k)])), (self.p, 1))
        self.y_tau = 1.0 / (1.0 - h2)  # var_Y get cancelled eventually
        self.prior_pi = np.ones((self.p,)) * (1 / self.p)
        self.beta_post_tau = np.tile(XX.reshape(-1, 1), (1, self.k)) * self.y_tau + self.beta_prior_tau

    def infer_q_beta(self, ytX, XtX):
        """perform variational updates"""
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            beta_all_k = (self.gamma[:, idxall] * self.beta_mu[:, idxall]).sum(axis=1)
            self.beta_mu[:, k] = (ytX - np.dot(beta_all_k, XtX)) / self.beta_post_tau[:, k] * self.y_tau  # update mu
            u = -0.5 * np.log(self.beta_post_tau[:, k]) + np.log(self.prior_pi.transpose()) + \
                0.5 * self.beta_mu[:, k] ** 2 * self.beta_post_tau[:, k]
            self.gamma[:, k] = softmax(u)  # update gamma

    def get_elbo(self, XX, ytX, XtX):
        beta_all = (self.gamma * self.beta_mu).sum(axis=1)
        ll1 = - 2 * np.dot(beta_all, ytX)
        ll2 = ((self.gamma * self.beta_mu ** 2).sum(axis=1) * XX).sum()
        W = self.gamma * self.beta_mu
        WtRW = np.dot(np.dot(W.transpose(), XtX), W)
        ll3 = WtRW.sum() - np.diag(WtRW).sum()
        ll = -0.5 * self.y_tau * (ll1 + ll2 + ll3)
        betaterm1 = -0.5 * (self.beta_prior_tau * self.gamma * (self.beta_mu ** 2)).sum()
        gammaterm1 = (self.gamma * np.tile(self.prior_pi.reshape(-1, 1), (1, self.k))).sum()
        gammaterm2 = (self.gamma[self.gamma != 0] * np.log(self.gamma[self.gamma != 0])).sum()
        eut = -0.5 * (self.gamma * np.log(self.beta_post_tau)).sum()  # extra unstandardized term
        mkl = betaterm1 + gammaterm1 - gammaterm2 + eut
        elbo = ll + mkl
        return ll, mkl, elbo

    def get_PIP(self):
        return np.max((self.gamma), axis=1)

    def update_pi(self, new_pi):
        self.prior_pi = new_pi

    def get_effect(self, cthres=0.95, ethres=20):
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k))  # effective gamma
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        csum = mat_eff.sum(axis=0).round(2)
        print("Attainable coverage for effect groups: {}".format(csum))  # statistical evidence
        eff = {}
        eff_gamma = {}
        eff_mu = {}
        for k in range(self.k):
            if csum[k] > cthres:
                if entropy(mat_eff[:, k]) < np.log(ethres):
                    for p in range(self.p):
                        if np.sum(mat_eff[matidx[0:p, k], k]) > cthres * csum[k] or mat_eff[matidx[p, k], k] < 0.01:
                            eff[k] = matidx[0:p, k]
                            eff_gamma[k] = mat_eff[eff[k], k].round(4)
                            eff_mu[k] = self.beta_mu[eff[k], k].round(4)
                            break
        return eff, eff_gamma, eff_mu

    def train(self, XX, ytX, XtX, maxite=100, eps=0.01, verbose=False, loss=0.0):
        for ite in range(maxite):
            self.infer_q_beta(ytX, XtX)
            ll, mkl, elbo = self.get_elbo(XX, ytX, XtX)
            if verbose:
                print('*' * 70)
                print('Iteration-->{} . Likelihood: {:.2f} . KL: {:.2f} . ELBO: {:.2f}'.format(ite, ll, mkl, elbo))
            if abs(elbo - loss) < eps:
                break
            if ite == (maxite - 1):
                print("Algorithm not converged. Please make sure matched summary statistics and LD were provided!")
            loss = elbo


def get_XX_XtX_ytX_z(LD, Z, N):
    assert not np.isnan(LD).any(), "The LD matrix contains NaN values."
    assert not np.isnan(Z).any(), "The zscore files contains NaN values."
    XX = np.ones(len(Z)) * N
    XtX = LD * N
    ytX = Z * np.sqrt(N)
    return XX, XtX, ytX


def get_HESS_h2_z(LD, Z, N, ptLD=0.2, ptp=1e-5):
    """calculate local heritabilities"""
    assert not np.isnan(LD).any(), "The LD matrix contains NaN values."
    assert not np.isnan(Z).any(), "The zscore files contains NaN values."
    zsquare = Z ** 2
    pvalmin = chi2.sf(max(zsquare), 1)
    idx_retain = []
    idx_exclude = [i for i in range(len(Z))]
    while len(idx_exclude) > 0:  # P + T
        maxid = idx_exclude[np.argmax(zsquare[idx_exclude])]  # find the idx with smallest p-value
        pval = chi2.sf(zsquare[maxid], 1)
        idx_retain.append(maxid)
        idx_exclude = [i for i in idx_exclude if i not in np.where(np.abs(LD[maxid, :]) > ptLD)[0]]
        if pval > ptp:
            break
    Indidx = np.sort(idx_retain)  # obtain independent signals
    P = len(Indidx)
    LD_id = LD[np.ix_(Indidx, Indidx)]
    R_inv = np.linalg.inv(LD_id)
    vec_id = Z[Indidx]
    h2_hess = (np.dot(np.dot(vec_id.transpose(), R_inv), vec_id) - P) / (N - P)
    var_b = np.max(zsquare[Indidx]) / N
    if h2_hess < 0.0001:
        h2_hess = 0.0001
    if h2_hess > 0.9:
        h2_hess = 0.9
    return h2_hess, var_b, pvalmin, P


def run_sparsepro(LD, cthres, ethres, Z, N, K, h2_hess, var_b, verbose, new_pi=None):
    XX, XtX, ytX = get_XX_XtX_ytX_z(LD, Z, N)
    model = SparsePro(len(Z), K, XX, h2_hess, var_b)
    if new_pi is not None:
        model.update_pi(new_pi)
    model.train(XX, ytX, XtX, verbose=verbose)
    mcs, eff_gamma, eff_mu = model.get_effect(cthres, ethres)
    pip_vec = model.get_PIP().round(4)
    return mcs, eff_gamma, eff_mu, pip_vec


def chi_squared(row):
    if any(row[0:4] == 0):
        return '1.00e+00'
    chi2, p, dof, expctd = chi2_contingency(np.array([[row[0], row[1]], [row[2], row[3]]]), lambda_="log-likelihood")
    return '{:.2e}'.format(p)


def get_uni_enrich(annomat, pip):
    df = pd.DataFrame(0.0, columns=['k0', 'k1', 'r0', 'r1'], index=annomat.columns)
    for a in annomat.columns:
        df.loc[a, 'k0'] = (1 - annomat[a]).sum()
        df.loc[a, 'k1'] = annomat[a].sum()
        df.loc[a, 'r0'] = np.dot(pip, (1 - annomat[a]).values)
        df.loc[a, 'r1'] = np.dot(pip, annomat[a].values)
    df['W'] = np.log(df['r1'] * df['k0'] / df['r0'] / df['k1'])
    df['W_se'] = np.sqrt(1 / df['r1'] + 1 / df['r0'] - 1 / df['k1'] - 1 / df['k0'])
    df = df.round(2)
    df['p'] = df.apply(chi_squared, axis=1)
    return df


def get_all_enrich(annomat, pip, maxite=20):
    W_new = pd.Series(0.0, index=annomat.columns)
    for ite in range(maxite):
        W_old = W_new.copy()
        df = pd.DataFrame(0.0, columns=['k0', 'k1', 'r0', 'r1'], index=annomat.columns)
        for a in annomat.columns:
            wt = softmax(np.dot(annomat.drop(a, axis=1).values, W_new.drop(a))).squeeze()
            df.loc[a, 'k0'] = np.dot(wt, (1 - annomat[a]).values) * annomat.shape[0]
            df.loc[a, 'k1'] = np.dot(wt, annomat[a].values) * annomat.shape[0]
            df.loc[a, 'r0'] = np.dot(pip, (1 - annomat[a]).values)
            df.loc[a, 'r1'] = np.dot(pip, annomat[a].values)
            W_new[a] = np.log(df.loc[a, 'r1'] * df.loc[a, 'k0'] / df.loc[a, 'r0'] / df.loc[a, 'k1'])
        eps = ((W_new - W_old) ** 2).sum()
        print("iteration {} with diff {}".format(ite, eps))
        if eps < 1e-2:
            break
    df = df.round(2)
    df['W'] = W_new.round(2)
    df['W_se'] = np.sqrt(1 / df['r1'] + 1 / df['r0'] - 1 / df['k1'] - 1 / df['k0']).round(2)
    df['p'] = df.apply(chi_squared, axis=1)
    return df


def zld(args):
    pip = []
    ldlists = pd.read_csv(args.zld, sep='\s+')  # z ld anno
    print("LD list with {} LD blocks loaded\n".format(len(ldlists)))
    if args.aW:
        W_sig = pd.read_csv(args.aW, sep='\t')
        sigidx = W_sig['sigidx'].tolist()
        W_new = W_sig['W'].values
    for ite in range(len(ldlists)):
        cs = []
        cs_pip = []
        cs_eff = []
        zfile, ldfile = ldlists.iloc[ite, 0:2]
        LD = pd.read_csv(os.path.join(args.zdir, ldfile), sep='\s+', header=None)  # no header square symmetric matrix
        z = pd.read_csv(os.path.join(args.zdir, zfile), sep='\t', header=None, index_col=0)  # no header matched vector
        assert len(LD) == len(z), 'LD file and Z file must match!'
        if args.h2:
            h2_hess, var_b, p_chi, K = ldlists.loc[ite, ['h2', 'varb', 'pval', 'K']]
        else:
            h2_hess, var_b, p_chi, K = get_HESS_h2_z(LD.values, z[1].values, args.N)
        if args.K is not None:
            K = args.K
        if args.anno is not None and args.aW is not None:
            annofile = ldlists.at[ite, args.anno]
            ianno = pd.read_csv(os.path.join(args.zdir, annofile), sep='\t', index_col=0)
            ANN = ianno.values[:, sigidx]
            new_pi = softmax(np.dot(ANN, W_new))
            mcs, eff_gamma, eff_mu, pip_vec = run_sparsepro(LD.values, args.cthres, args.ethres,
                                                            z[1].values, args.N, K, h2_hess,
                                                            var_b, args.verbose, new_pi)
        else:
            mcs, eff_gamma, eff_mu, pip_vec = run_sparsepro(LD.values, args.cthres, args.ethres,
                                                            z[1].values, args.N, K, h2_hess,
                                                            var_b, args.verbose)
        print("{} variants loaded from {} explaining {:2.2%} of trait heritability with smallest p-value {:.2e}".
              format(LD.shape[1], zfile, h2_hess, p_chi))
        pip.extend(pip_vec)
        z['pip'] = pip_vec
        if len(mcs) == 0:
            print("No effect detected")
        else:
            print("Detected k = {}\n".format(len(mcs)))
            for e in mcs:
                mcs_idx = [z.index[j] for j in mcs[e]]
                print('The {}-th effect contains effective variants:'.format(e))
                print('causal variants: {}'.format(mcs_idx))
                print('posterior inclusion probabilities: {}'.format(eff_gamma[e]))
                print('posterior causal effect size: {}'.format(eff_mu[e]))
                print()
                cs.append(mcs_idx)
                cs_pip.append(eff_gamma[e])
                cs_eff.append(eff_mu[e])
        if not args.h2:
            ldlists.at[ite, 'h2'] = '{:.2e}'.format(h2_hess)
            ldlists.at[ite, 'pval'] = '{:.2e}'.format(p_chi)
            ldlists.at[ite, 'varb'] = '{:.2e}'.format(var_b)
            ldlists.at[ite, 'K'] = '{:.2e}'.format(K)
        z.to_csv(os.path.join(args.save, "{}.pip".format(zfile)), sep='\t', header=False, index=True)
        allcs = pd.DataFrame({"cs": ['/'.join(i) for i in cs], "pip": ['/'.join([str(j) for j in i]) for i in cs_pip],
                              "beta": ['/'.join([str(j) for j in i]) for i in cs_eff]})
        allcs.to_csv(os.path.join(args.save, "{}.cs".format(zfile)), sep='\t', header=True, index=False)
    ldlists.to_csv(os.path.join(args.save, "{}.h2".format(args.prefix)), sep='\t', header=True, index=False)

    if args.anno is not None and args.aW is None:
        annomat = pd.concat([pd.read_csv(os.path.join(args.zdir, i), sep='\s+', index_col=0)
                             for i in ldlists[args.anno]])
        df = get_uni_enrich(annomat, pip)
        sigidx = [i for i in range(annomat.shape[1]) if float(df['p'][i]) < args.pthres]
        df.reset_index().to_csv(os.path.join(args.save, '{}.wsep'.format(args.prefix)),
                                sep="\t", index=False)
        valididx = [i for i in range(annomat.shape[1]) if annomat.sum()[i]>0]
        validANNOT = annomat[annomat.columns[valididx]]
        df_all = get_all_enrich(validANNOT, pip)
        df_all['sigidx'] = valididx
        df_all.reset_index().to_csv(os.path.join(args.save, '{}.W1.0'.format(args.prefix)),
                                    sep='\t', index=False)
        if len(sigidx) == 0:
            print("None of the {} annotations is significantly enriched at p-value threshold {}.".
                  format(annomat.shape[1], args.pthres))
        else:
            print("{} annotations are deemed significantly enriched at {} p-value threshold and used to update priors."
                  "Saving result to {}.W{} file...".format(len(sigidx), args.pthres, args.prefix, args.pthres))
            sigANNOT = annomat[annomat.columns[sigidx]]
            df = get_all_enrich(sigANNOT, pip)
            df['sigidx'] = sigidx
            df.reset_index().to_csv(os.path.join(args.save, '{}.W{}'.format(args.prefix, args.pthres)),
                                    sep='\t', index=False)


parser = argparse.ArgumentParser(description='SparsePro Commands:')
parser.add_argument('--zld', type=str, default=None,
                    help='locus finemapping mode: path to matched zscore and ld lists (or annotation file)')
parser.add_argument('--zdir', type=str, default=None, help='path to zscores files', required=True)
parser.add_argument('--N', type=int, default=None, help='GWAS sample size', required=True)
parser.add_argument('--save', type=str, default=None, help='path to save result', required=True)
parser.add_argument('--prefix', type=str, default=None, help='prefix for result files', required=True)
parser.add_argument('--verbose', action="store_true", help='options for displaying more information')
parser.add_argument('--anno', type=str, default=None, help='name of the annotation column in zld file')
parser.add_argument('--K', type=int, default=None, help='largest number of effect groups')
parser.add_argument('--pthres', type=float, default=1e-5, help='p value threshold for enrichment')
parser.add_argument('--cthres', type=float, default=0.95, help='coverage level for credible sets')
parser.add_argument('--ethres', type=float, default=20.0, help='entropy level for credible sets')
parser.add_argument('--aW', type=str, default=None, help='significant enriched file')
parser.add_argument('--h2', action="store_true", help='use previous h2 file as zld file')

args = parser.parse_args()
title()
if not os.path.exists(args.save):
    os.makedirs(args.save)
zld(args)
