import numpy as np

def test_preprocess_random_data(tmp_path):

    import loompy

    from panopticon.preprocessing import generate_count_normalization
    from panopticon.analysis import generate_incremental_pca, generate_embedding, generate_clustering, generate_masked_module_score

    d = tmp_path / "sub"
    d.mkdir()
    p = str(d / "test.loom")
    data = np.random.randint(1000,size=(100,1000))
    loompy.create(p,data,row_attrs = {'gene': np.arange(100)},col_attrs = {'cell':np.arange(1000)})
    db = loompy.connect(p)

    generate_count_normalization(db, '','log2(TP100k+1)')
    generate_incremental_pca(db, 'log2(TP100k+1)')
    generate_embedding(db, 'log2(TP100k+1)' )
    generate_clustering(db, 'log2(TP100k+1)', clusteringcachedir=str(d/"clusteringcachedir/"))
    print(db.ca.keys())
    print(db.ra.keys())
    assert 'cell' in db.ca.keys()
    assert 'gene' in db.ra.keys()
    assert 'PCA UMAP embedding 1' in db.ca.keys()
    db.close()



