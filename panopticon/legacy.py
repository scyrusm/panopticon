import os

import numpy as np
import pandas as pd


def get_module_score_loom(loom,
                          signature_name,
                          querymask=None,
                          nbins=10,
                          ncontrol=5):
    """Calculates a module score over a loom file.  This routine is deprecated--use generate masked module score (S Markson 3 June 2020).

    .. deprecated:: 0.1 

    Parameters
    ----------
    loom : loom object on which to calculate score
        
    signature_name : Name of signature (pre-loaded into loom object) over which to calculate score
        
    nbins : Number of quantile bins to use
        (Default value = 100)
    ncontrol : Number of genes in each matched quantile
        (Default value = 5)
    querymask :
        (Default value = None)

    Returns
    -------

    
    """
    if querymask is None:
        querymask = np.array([True] * loom.shape[1])
    if querymask.dtype == bool:
        querymask = querymask.nonzero()[0]

    alldata = loom[:, querymask]
    nonsigdata = alldata[loom.ra[signature_name] == 0]
    sigdata = alldata[loom.ra[signature_name] == 1]

    gene_quantiles = pd.qcut(alldata.mean(axis=1),
                             nbins,
                             duplicates='drop',
                             labels=False)
    sigdata_quantiles = gene_quantiles[loom.ra[signature_name] == 1]
    nonsigdata_quantiles = gene_quantiles[loom.ra[signature_name] == 0]
    signature = sigdata.mean(axis=0)
    control_group = []
    for quantile in np.unique(sigdata_quantiles):
        noccurrences = np.sum(sigdata_quantiles == quantile)
        # there's an edge case wherein if size is greater than the number of genes to be taken without replacement, this will generate an error.  Will be an issue in the few-gene regime
        control_group += list(
            np.random.choice(np.where(nonsigdata_quantiles == quantile)[0],
                             size=ncontrol * noccurrences,
                             replace=False))

    control_group = np.array(control_group)
    control = nonsigdata[control_group].mean(axis=0)
    return signature - control


def create_subsetted_loom_with_genemask(loom, output_loom, cellmask, genemask):
    """Deprecated.

    Parameters
    ----------
    loom :
        
    output_loom :
        
    cellmask :
        
    genemask :
        

    Returns
    -------

    
    """
    print("THIS FUNCTION IS DEPRECATED, USE loompy.new INSTEAD!!!")
    import loompy

    from panopticon.utilities import recover_meta
    if '' not in loom.layers.keys():
        raise Exception("Expecting '' layer, yet none found")
    rowmeta, colmeta = recover_meta(loom)
    if len(genemask)!=loom.shape[0] or len(cellmask)!=loom.shape[1]:
        raise Exception("genemask and cellmask must be boolean masks with length equal to the number of rows and columns of loom, respectively")
    loompy.create(output_loom, loom[''][genemask.nonzero()[0], :][:, cellmask.nonzero()[0]],
                  rowmeta[genemask].to_dict("list"),
                  colmeta[cellmask].to_dict("list"))
    with loompy.connect(output_loom) as smallerloom:
        for layer in [x for x in loom.layer.keys() if x != '']:
            smallerloom[layer] = loom[layer][:, cellmask.nonzero()[0]][genemask.nonzero()[0], :]


def create_subsetted_loom(loom, output_loom_filename, cellmask):
    """Deprecated.

    Will create a new loom file with cells specified according to a Boolean vector mask.

    Parameters
    ----------
    loom : LoomConnection object which will be subsetted
        
    output_loom_filename : string denoting the path and filename of the output loom file.  
        
    cellmask : Boolean numpy vector with length equal to the number of cells in "loom"
        
    Returns
    -------

    
    """
    import loompy

    if len(cellmask)!=loom.shape[1]:
        raise Exception("cellmask must be boolean mask with length equal to the number of columns of loom")
    with loompy.new(output_loom_filename) as dsout:
        cells = np.where(cellmask)[
            0]  
        for (ix, selection, view) in loom.scan(items=cells, axis=1,
                                               key="gene"):
            dsout.add_columns(view.layers,
                              col_attrs=view.ca,
                              row_attrs=view.ra)


def get_gsea_with_selenium(diffex,
                           email='s'):
    """If you aren't Sam, probably don't use this.

    Parameters
    ----------
    diffi :
        

    Returns
    -------

    
    """

    from selenium import webdriver
    from selenium.webdriver.support.ui import Select
    profile = webdriver.FirefoxProfile()
    profile.set_preference('browser.download.folderList', 2) # custom location
    profile.set_preference('browser.download.manager.showWhenStarting', False)
    profile.set_preference('browser.download.dir', '/tmp')
    profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/tsv')
    profile.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/txt')
    
    driver = webdriver.Firefox(profile)
    download_folder = '/tmp/'
    genelisthitdict = {}
    flag = True
    email = 's'
    for key in diffex.keys():
    
        driver.get("https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp")
        if flag:
            elem = driver.find_element_by_name("j_username")
            elem.send_keys(email)
            elem = driver.find_element_by_name('login')
            elem.click()
            flag = False
        elem = driver.find_element_by_id("geneList")
        elem.send_keys('\n'.join(diffex[key].query('MeanExpr1 > MeanExpr2')['gene'][0:200].values))
        for genelistcheckbox in ["H", "C2", "C5", "C6", "C7","C8"]:
            elem = driver.find_element_by_id(
                "{}checkbox".format(genelistcheckbox))
            elem.click()
        #driver.close()
        #elem.find_element_by_xpath("/html/body/div[2]/div[4]/div/table/tbody/tr/td[3]/table/tbody/tr/td[2]/form/table/tbody/tr[30]/td/p[1]/nobr/select[@name='element_name']/option[text()='100']")
        #elem.click()
        select = Select(driver.find_element_by_xpath('/html/body/div[2]/div[4]/div/table/tbody/tr/td[3]/table/tbody/tr/td[2]/form/table/tbody/tr[30]/td/p[1]/nobr/select'))
        select.select_by_value('100')
    
        elem = driver.find_element_by_name("viewOverlaps")
        elem.click()
        #try:
        elem = driver.find_element_by_link_text("Text")
        elem.click()
        #driver.close()
        footersize = 0
        while True:
            try:
                genelisthits = pd.read_table(
                    "{}/overlap.tsv".format(download_folder),
                    skiprows=10,
                    skipfooter=footersize,
                    header=None)
                break
            except:
                footersize += 1
        genelisthits = pd.read_table(
            "{}/overlap.tsv".format(download_folder),
            skiprows=10,
            skipfooter=footersize + 2,
            header=None)
        genelisthitdict[key] = genelisthits
        os.system("rm {}/overlap.tsv".format(download_folder))
        #except:
        #    genelisthitdict[key] = 'no hits'
        #print(genelisthitdict[key])
    driver.close()
    return genelisthitdict
