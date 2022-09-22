import numpy as np
import pandas as pd


def gravy_index(cdr3):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return [ProteinAnalysis(x.replace('-', '')).gravy() for x in cdr3]


def multinomial_test(x1, x2, x12, n, doublet_rate):
    c = doublet_rate / (1 - doublet_rate) / (1 - doublet_rate)
    from scipy.special import comb

    def g(x1, x2, x12, n, c):
        return c**x12 * comb(x1 + x12, x12, exact=False) * comb(
            n - x2 - x12, x1, exact=False)

    numerator = np.sum([
        g(x1 + x12 - x, x2 + x12 - x, x, n, c)
        for x in range(x12, 1 + np.min([x1 + x12, x2 + x12]))
    ])
    denominator = np.sum([
        g(x1 + x12 - x, x2 + x12 - x, x, n, c)
        for x in range(1, 1 + np.min([x1 + x12, x2 + x12]))
    ])
    return numerator / denominator


def get_tcr_summary_df(tcrs,
                       samples=None,
                       multichain_separator='---',
                       alpha_beta_separator='|',
                       sample_separator='_',
                       tcrblacklist=[],
                       multinomial_test=True):

    from panopticon.tcr import get_tcr_doublet_likelihood

    df = pd.DataFrame(tcrs, copy=True, columns=['TCR'])
    if samples is None:
        df['sample'] = 'full dataset'
    else:
        df['sample'] = samples
    df = df.groupby('sample')['TCR'].value_counts()
    df.name = 'counts'
    df = df.reset_index()
    #    df['frequencies'] = df.groupby('sample')['TCR'].value_counts(
    #        normalize=True).values

    df = df[~df['TCR'].isin(tcrblacklist)]
    df['TRA'] = df['TCR'].apply(
        lambda x: x.split(alpha_beta_separator)[0]
    )  # This implies certain expectations about how the TCR string is written
    df['TRB'] = df['TCR'].apply(
        lambda x: x.split(alpha_beta_separator)[1].split(sample_separator)[0]
    )  # string after the "_" indicates the sample, or another annotation

    df['N_TRA'] = df['TRA'].apply(
        lambda x: 0 if x == '' else len(x.split(multichain_separator)))
    df['N_TRB'] = df['TRB'].apply(
        lambda x: 0 if x == '' else len(x.split(multichain_separator)))
    if multinomial_test is False:
        return df

    a1b2 = get_tcr_doublet_likelihood(
        tcrs,
        samples=samples,
        tcrblacklist=tcrblacklist,
        multichain_separator=multichain_separator,
        alpha_beta_separator=alpha_beta_separator,
        sample_separator=sample_separator,
        n_tra=1,
        n_trb=2)
    a2b1 = get_tcr_doublet_likelihood(
        tcrs,
        samples=samples,
        tcrblacklist=tcrblacklist,
        multichain_separator=multichain_separator,
        alpha_beta_separator=alpha_beta_separator,
        sample_separator=sample_separator,
        n_tra=2,
        n_trb=1)
    a2b2 = get_tcr_doublet_likelihood(
        tcrs,
        samples=samples,
        tcrblacklist=tcrblacklist,
        multichain_separator=multichain_separator,
        alpha_beta_separator=alpha_beta_separator,
        sample_separator=sample_separator,
        n_tra=2,
        n_trb=2)
    aMbN = pd.concat([a1b2, a2b1, a2b2])
    df = pd.merge(df,
                  aMbN[[x for x in aMbN.columns if x not in df.columns] +
                       ['sample', 'TCR']],
                  on=['sample', 'TCR'],
                  how='outer',
                  validate='one_to_one')

    return df


def get_tcr_doublet_likelihood(tcrs,
                               samples=None,
                               tcrblacklist=[],
                               n_tra=2,
                               n_trb=2,
                               verbose=False,
                               multichain_separator='---',
                               alpha_beta_separator='|',
                               sample_separator='_',
                               doublet_rate=0.05):
    from panopticon.tcr import multinomial_test
    from panopticon.tcr import get_tcr_summary_df
    from statsmodels.stats.multitest import fdrcorrection
    from tqdm import tqdm
    df = get_tcr_summary_df(tcrs,
                            samples=samples,
                            tcrblacklist=tcrblacklist,
                            multichain_separator=multichain_separator,
                            alpha_beta_separator=alpha_beta_separator,
                            sample_separator=sample_separator,
                            multinomial_test=False)
    ps = []
    #n =
    pairs = []
    aMbN = df.query('N_TRA==@n_tra').query('N_TRB==@n_trb')
    for i, row in tqdm(aMbN.iterrows(),
                       desc='Looping over {} TCR-A {} TCR-B clonotypes'.format(
                           n_tra, n_trb),
                       total=len(aMbN)):
        sample = row['sample']
        putative_ps = []
        putative_pairs = []
        if n_tra == 2 and n_trb == 2:
            for alpha in row['TRA'].split(multichain_separator):
                alpha_c = [
                    x for x in row['TRA'].split(multichain_separator)
                    if x != alpha
                ][0]
                for beta in row['TRB'].split(
                        multichain_separator)[0:1]:  # Only want the first beta
                    beta_c = [
                        x for x in row['TRB'].split(multichain_separator)
                        if x != beta
                    ][0]
                    df1 = df.query('sample==@sample').query(
                        'TRA==@alpha').query('TRB==@beta')
                    df2 = df.query('sample==@sample').query(
                        'TRA==@alpha_c').query('TRB==@beta_c')
                    if df1.shape[0] == 1:
                        x1 = df1['counts'].values[0]
                    elif df1.shape[0] > 1:
                        raise Exception("!")
                    else:
                        x1 = 0
                    if df2.shape[0] == 1:
                        x2 = df2['counts'].values[0]
                    elif df2.shape[0] > 1:
                        raise Exception("!")
                    else:
                        x2 = 0
                    x12 = row['counts']

                    n = df.query('sample==@sample')['counts'].sum()
                    if x1 > 0 and x2 > 0:
                        p = multinomial_test(x1, x2, x12, n, doublet_rate)
                    else:
                        p = np.nan

                    putative_ps.append(p)
                    pair = '{}|{} (n={}) + {}|{} (n={})'.format(
                        alpha, beta, x1, alpha_c, beta_c, x2)
                    putative_pairs.append(pair)
        if n_tra == 2 and n_trb == 1:
            alpha, alpha_c = row['TRA'].split(multichain_separator)
            for beta in [row['TRB'], '*']:
                if beta == '*':
                    beta_c = row['TRB']
                    df1 = df.query('sample==@sample').query(
                        'TRA==@alpha').query('N_TRA==1').query('N_TRB==1')
                    df2 = df.query('sample==@sample').query(
                        'TRA==@alpha_c').query('TRB==@beta').query(
                            'N_TRA==1').query('N_TRB==1')
                else:
                    beta_c = '*'
                    df1 = df.query('sample==@sample').query(
                        'TRA==@alpha').query('TRB==@beta').query(
                            'N_TRA==1').query('N_TRB==1').sort_values(
                                'counts', ascending=False)
                    df2 = df.query('sample==@sample').query(
                        'TRA==@alpha_c').query('N_TRA==1').query(
                            'N_TRB==1').sort_values('counts', ascending=False)
                if df1.shape[0] == 1:
                    x1 = df1['counts'].values[0]
                elif df1.shape[0] > 1:
                    if verbose:
                        display(df1)
                    x1 = df1['counts'].sum()
                    #raise Exception("Multiple clones matching condition?")
                else:
                    x1 = 0
                if df2.shape[0] == 1:
                    x2 = df2['counts'].values[0]
                elif df2.shape[0] > 1:
                    if verbose:
                        display(df2)
                    x2 = df2['counts'].sum()
                # raise Exception("Multiple clones matching condition?")
                else:
                    x2 = 0
                x12 = row['counts']

                n = df.query('sample==@sample')['counts'].sum()
                if x1 > 0 and x2 > 0:
                    p = multinomial_test(x1, x2, x12, n, doublet_rate)
                else:
                    p = np.nan

                putative_ps.append(p)
                pair = '{}|{} (n={}) + {}|{} (n={})'.format(
                    alpha, beta, x1, alpha_c, beta_c, x2)
                putative_pairs.append(pair)
        if n_tra == 1 and n_trb == 2:
            beta, beta_c = row['TRB'].split(multichain_separator)
            for alpha in [row['TRA'], '*']:
                if alpha == '*':
                    alpha_c = row['TRA']
                    df1 = df.query('sample==@sample').query(
                        'TRA==@alpha').query('TRB==@beta').query(
                            'N_TRA==1').query('N_TRB==1').sort_values(
                                'counts', ascending=False)
                    df2 = df.query('sample==@sample').query(
                        'TRB==@beta_c').query('N_TRA==1').query(
                            'N_TRB==1').sort_values('counts', ascending=False)
                else:
                    alpha_c = '*'
                    df1 = df.query('sample==@sample').query(
                        'TRB==@beta').query('N_TRA==1').query(
                            'N_TRB==1').sort_values('counts', ascending=False)
                    df2 = df.query('sample==@sample').query(
                        'TRA==@alpha').query('TRB==@beta_c').query(
                            'N_TRA==1').query('N_TRB==1').sort_values(
                                'counts', ascending=False)
                if df1.shape[0] == 1:
                    x1 = df1['counts'].values[0]
                elif df1.shape[0] > 1:
                    if verbose:
                        display(df1)
                    x1 = df1['counts'].sum()
                #  raise Exception("Multiple clones matching condition?")
                else:
                    x1 = 0
                if df2.shape[0] == 1:
                    x2 = df2['counts'].values[0]
                elif df2.shape[0] > 1:
                    if verbose:
                        display(df2)
                    x2 = df2['counts'].sum()
                #    raise Exception("Multiple clones matching condition?")
                else:
                    x2 = 0
                x12 = row['counts']

                n = df.query('sample==@sample')['counts'].sum()
                if x1 > 0 and x2 > 0:
                    p = multinomial_test(x1, x2, x12, n, doublet_rate)
                else:
                    p = np.nan

                putative_ps.append(p)
                pair = '{}|{} (n={}) + {}|{} (n={})'.format(
                    alpha, beta, x1, alpha_c, beta_c, x2)
                putative_pairs.append(pair)
        if np.sum(~np.isnan(putative_ps)) == 0:
            ps.append(np.nan)
            pairs.append('No putative pairs with positive population')
        else:
            arg = np.nanargmax(putative_ps)
            ps.append(putative_ps[arg])
            pairs.append(putative_pairs[arg])
    aMbN['Multinomial_p'] = ps
    aMbN['MostLikelyDoublePair'] = pairs
    nonnanp = [x for x in ps if not np.isnan(x)]
    relevantq = fdrcorrection(nonnanp, is_sorted=False)[1]
    p2q = {p: q for p, q in zip(nonnanp, relevantq)}
    aMbN['Multinomial_BenjaminiHochbergQ'] = [
        p2q[p] if not np.isnan(p) else np.nan
        for p in aMbN['Multinomial_p'].values
    ]
    return aMbN


def generate_tcr_multichain_summary(loom,
                                    tcr_ca='TCR',
                                    sample_ca=None,
                                    tcrblacklist=[],
                                    suffix='_multichain_summary',
                                    multichain_separator='---',
                                    alpha_beta_separator='|',
                                    sample_separator='_',
                                    overwrite=False):
    from panopticon.tcr import get_tcr_summary_df
    if sample_ca is None:
        samples = None
    else:
        samples = loom.ca[sample_ca]
    df = get_tcr_summary_df(loom.ca[tcr_ca],
                            samples=samples,
                            multichain_separator=multichain_separator,
                            alpha_beta_separator=alpha_beta_separator,
                            sample_separator=sample_separator,
                            tcrblacklist=tcrblacklist)
    dfdict = df.set_index(['sample', 'TCR']).to_dict()
    for key in dfdict.keys():
        if key + suffix in loom.ca.keys() and overwrite is False:
            raise Exception(
                "{} already in column attribute keys of {}; adjust suffix or set overwrite to True"
                .format(key + suffix, loom.filename))
        else:
            if sample_ca is None:
                iterator = zip(['full dataset'] * loom.shape[1],
                               loom.ca[tcr_ca])
            else:
                iterator = zip(loom.ca[sample_ca], loom.ca[tcr_ca])
            loom.ca[key + suffix] = [
                dfdict[key][(sample, tcr)] if
                (sample, tcr) in dfdict[key].keys() else np.nan
                for sample, tcr in iterator
            ]


def incorporate_10x_vdj(loomfile,
                        filtered_contig_annotations_csv,
                        barcode_ca='cellname',
                        overwrite=False,
                        barcode_match_exception_threshold=0.5):
    """

    Parameters
    ----------
    loom :
        
    filtered_contig_annotations_csv :
        
    barcode_ca :
         (Default value = 'cellname')

    Returns
    -------

    """
    import loompy
    loom = loompy.connect(loomfile)
    if len(np.unique(loom.ca[barcode_ca])) != loom.shape[1]:
        raise Exception("`barcode_ca` must be unique to each cell")

    filtered_contig_annotations = pd.read_csv(filtered_contig_annotations_csv)
    barcode_match_rate = np.isin(
        filtered_contig_annotations['barcode'].unique(),
        loom.ca[barcode_ca]).mean()
    if barcode_match_rate < barcode_match_exception_threshold:
        raise Exception(
            "Only {}% of V(D)J barcodes have corresponding gene expression data; GEX and V(D)J files may be mismatched"
            .format(100 * barcode_match_rate))
    filtered_contig_annotations_tcra = filtered_contig_annotations.query(
        'chain=="TRA"').query('productive=="True"')
    filtered_contig_annotations_tcrb = filtered_contig_annotations.query(
        'chain=="TRB"').query('productive=="True"')
    barcode_df_dict = {}
    for label, fca in zip(
        ['TRA', 'TRB'],
        [filtered_contig_annotations_tcra, filtered_contig_annotations_tcrb]):
        barcode_df = None
        for col in [
                x for x, xtype in zip(fca.columns, fca.dtypes)
                if x != 'barcode' and xtype not in [bool, int, float]
        ]:
            if (~fca[col].isnull()).sum() > 0:
                fca[col] = fca[col].astype(str)
                if barcode_df is None:
                    barcode_df = fca.groupby('barcode').agg({
                        col: '---'.join
                    }).reset_index()
                else:
                    barcode_df = pd.merge(barcode_df,
                                          fca.groupby('barcode').agg({
                                              col:
                                              '---'.join
                                          }).reset_index(),
                                          on='barcode')
        barcode_df_dict[label] = barcode_df.set_index('barcode').to_dict()
        #break
    for superkey in barcode_df_dict:
        for subkey in barcode_df_dict[superkey]:
            if '_'.join([superkey, subkey
                         ]) in loom.ca.keys() and overwrite == False:
                raise Exception(
                    '{} is already a column attribute in {}; cannot re-assign while `overwrite` is False'
                    .format('_'.join([superkey, subkey]), loom.filename))
            else:
                loom.ca['_'.join([superkey, subkey])] = [
                    barcode_df_dict[superkey][subkey][x] if x
                    in barcode_df_dict[superkey][subkey].keys() else np.nan
                    for x in loom.ca[barcode_ca]
                ]


def join_tra_trb_ca(loom, ca='cdr3'):
    loom.ca['TRA_TRB_{}'.format(ca)] = [
        '|'.join([x, y]) for x, y, in zip(loom.ca['TRA_{}'.format(ca)],
                                          loom.ca['TRB_{}'.format(ca)])
    ]


def morisita(df, key, samplekey, sample1, sample2):
    from panopticon.analysis import simpson
    set1 = df[df[samplekey] == sample1][key].value_counts(
        normalize=True).reset_index(name='count', ).rename({'index': key},
                                                           axis=1)
    set2 = df[df[samplekey] == sample2][key].value_counts(
        normalize=True).reset_index(name='count', ).rename({'index': key},
                                                           axis=1)
    simpson1 = simpson(set1['count'].values, with_replacement=True)
    simpson2 = simpson(set2['count'].values, with_replacement=True)
    #allrearrangments = np.unique(np.hstack((set1['rearrangement'].values, set2['rearrangement'].values)))
    mergeset = pd.merge(set1,
                        set2,
                        on=key,
                        how='outer',
                        suffixes=('_set1', '_set2'))
    cross = (mergeset['count_set1'].fillna(0) *
             mergeset['count_set2'].fillna(0)).sum()
    return 2 * cross / (simpson1 + simpson2)


def _morisita(counts1, counts2):
    from panopticon.analysis import simpson
    simpson1 = simpson(counts1, with_replacement=True)
    simpson2 = simpson(counts2, with_replacement=True)
    cross = (counts1 * counts2).sum()/counts1.sum()/counts2.sum()
    return 2 * cross / (simpson1 + simpson2)
