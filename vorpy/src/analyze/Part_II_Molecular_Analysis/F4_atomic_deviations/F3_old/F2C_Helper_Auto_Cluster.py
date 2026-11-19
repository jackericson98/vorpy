import numpy as np
import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler



try:
    import hdbscan
    HAS_HDBSCAN = True

except ImportError:
    HAS_HDBSCAN = False


try:
    from kmodes.kprototypes import KPrototypes
    HAS_KPROTOTYPES = True

except ImportError:
    HAS_KPROTOTYPES = False



DEFAULT_NUMERIC_COLS = ['AW', 'Pow']
DEFAULT_CATEGORICAL_COLS = ['AtomName', 'ResidueName']
DEFAULT_BOOLEAN_COLS = []
DEFAULT_PERCENT_COLS = ['SolFacingPct']



def ensure_required_columns(
    df: pd.DataFrame,
    required_cols
):
    missing = [col for col in required_cols if col not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")



def make_sol_facing_binary(
    df: pd.DataFrame,
    source_col: str = 'SolFacingPct',
    out_col: str = 'SolFacingBinary',
    threshold: float = 20.0
) -> pd.DataFrame:
    df = df.copy()

    if source_col not in df.columns:
        raise ValueError(f"Column '{source_col}' not found in DataFrame.")

    df[out_col] = (df[source_col].astype(float) >= threshold).astype(int)

    return df



def build_feature_matrix(
    df: pd.DataFrame,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None,
    scale_numeric: bool = True,
    one_hot_categorical: bool = True
):
    if numeric_cols is None:
        numeric_cols = DEFAULT_NUMERIC_COLS

    if categorical_cols is None:
        categorical_cols = DEFAULT_CATEGORICAL_COLS

    if boolean_cols is None:
        boolean_cols = DEFAULT_BOOLEAN_COLS

    required_cols = list(numeric_cols) + list(categorical_cols) + list(boolean_cols)
    ensure_required_columns(df, required_cols)

    feature_blocks = []
    feature_names = []

    if numeric_cols:
        X_num = df[numeric_cols].astype(float).to_numpy()

        if scale_numeric:
            scaler = StandardScaler()
            X_num = scaler.fit_transform(X_num)

        feature_blocks.append(X_num)
        feature_names.extend(numeric_cols)

    if boolean_cols:
        X_bool = df[boolean_cols].astype(int).to_numpy()
        feature_blocks.append(X_bool)
        feature_names.extend(boolean_cols)

    if categorical_cols:
        X_cat_df = df[categorical_cols].astype(str)

        if one_hot_categorical:
            encoder = OneHotEncoder(
                sparse_output=False,
                handle_unknown='ignore'
            )

            X_cat = encoder.fit_transform(X_cat_df)
            cat_feature_names = encoder.get_feature_names_out(categorical_cols).tolist()

            feature_blocks.append(X_cat)
            feature_names.extend(cat_feature_names)

        else:
            for col in categorical_cols:
                codes = X_cat_df[col].astype('category').cat.codes.to_numpy().reshape(-1, 1)
                feature_blocks.append(codes)
                feature_names.append(col)

    if not feature_blocks:
        raise ValueError("No features were created. Check column inputs.")

    X = np.hstack(feature_blocks)

    return X, feature_names



def build_kprototypes_matrix(
    df: pd.DataFrame,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None,
    scale_numeric: bool = True
):
    if numeric_cols is None:
        numeric_cols = DEFAULT_NUMERIC_COLS

    if categorical_cols is None:
        categorical_cols = DEFAULT_CATEGORICAL_COLS

    if boolean_cols is None:
        boolean_cols = DEFAULT_BOOLEAN_COLS

    required_cols = list(numeric_cols) + list(categorical_cols) + list(boolean_cols)
    ensure_required_columns(df, required_cols)

    numeric_df = df[numeric_cols].astype(float).copy()

    if boolean_cols:
        for col in boolean_cols:
            numeric_df[col] = df[col].astype(int)

    X_num = numeric_df.to_numpy()

    if scale_numeric:
        scaler = StandardScaler()
        X_num = scaler.fit_transform(X_num)

    X_cat = df[categorical_cols].astype(str).to_numpy()

    X = np.hstack([X_num, X_cat])

    categorical_indices = list(range(X_num.shape[1], X_num.shape[1] + X_cat.shape[1]))

    return X, categorical_indices



def run_kmeans(
    df: pd.DataFrame,
    n_clusters: int = 8,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None,
    random_state: int = 0
) -> pd.DataFrame:
    out = df.copy()

    X, _ = build_feature_matrix(
        df=out,
        numeric_cols=numeric_cols,
        categorical_cols=categorical_cols,
        boolean_cols=boolean_cols,
        scale_numeric=True,
        one_hot_categorical=True
    )

    model = KMeans(
        n_clusters=n_clusters,
        n_init=20,
        random_state=random_state
    )

    out['MLCluster'] = model.fit_predict(X)

    return out



def run_agglomerative(
    df: pd.DataFrame,
    n_clusters: int = 8,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None,
    linkage: str = 'ward'
) -> pd.DataFrame:
    out = df.copy()

    X, _ = build_feature_matrix(
        df=out,
        numeric_cols=numeric_cols,
        categorical_cols=categorical_cols,
        boolean_cols=boolean_cols,
        scale_numeric=True,
        one_hot_categorical=True
    )

    model = AgglomerativeClustering(
        n_clusters=n_clusters,
        linkage=linkage
    )

    out['MLCluster'] = model.fit_predict(X)

    return out



def run_dbscan(
    df: pd.DataFrame,
    eps: float = 0.8,
    min_samples: int = 10,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None
) -> pd.DataFrame:
    out = df.copy()

    X, _ = build_feature_matrix(
        df=out,
        numeric_cols=numeric_cols,
        categorical_cols=categorical_cols,
        boolean_cols=boolean_cols,
        scale_numeric=True,
        one_hot_categorical=True
    )

    model = DBSCAN(
        eps=eps,
        min_samples=min_samples
    )

    out['MLCluster'] = model.fit_predict(X)

    return out



def run_hdbscan(
    df: pd.DataFrame,
    min_cluster_size: int = 25,
    min_samples: int = 10,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None
) -> pd.DataFrame:
    if not HAS_HDBSCAN:
        raise ImportError(
            "hdbscan is not installed. Install with: pip install hdbscan"
        )

    out = df.copy()

    X, _ = build_feature_matrix(
        df=out,
        numeric_cols=numeric_cols,
        categorical_cols=categorical_cols,
        boolean_cols=boolean_cols,
        scale_numeric=True,
        one_hot_categorical=True
    )

    model = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples
    )

    out['MLCluster'] = model.fit_predict(X)
    out['MLClusterProb'] = model.probabilities_

    return out



def run_kprototypes(
    df: pd.DataFrame,
    n_clusters: int = 8,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None,
    random_state: int = 0
) -> pd.DataFrame:
    if not HAS_KPROTOTYPES:
        raise ImportError(
            "kmodes is not installed. Install with: pip install kmodes"
        )

    out = df.copy()

    X, categorical_indices = build_kprototypes_matrix(
        df=out,
        numeric_cols=numeric_cols,
        categorical_cols=categorical_cols,
        boolean_cols=boolean_cols,
        scale_numeric=True
    )

    model = KPrototypes(
        n_clusters=n_clusters,
        init='Cao',
        n_init=5,
        random_state=random_state
    )

    out['MLCluster'] = model.fit_predict(
        X,
        categorical=categorical_indices
    )

    return out



def run_clustering(
    df: pd.DataFrame,
    method: str,
    numeric_cols=None,
    categorical_cols=None,
    boolean_cols=None,
    **kwargs
) -> pd.DataFrame:
    method = method.strip().lower()

    if method == 'kmeans':
        return run_kmeans(
            df=df,
            numeric_cols=numeric_cols,
            categorical_cols=categorical_cols,
            boolean_cols=boolean_cols,
            **kwargs
        )

    if method == 'agglomerative':
        return run_agglomerative(
            df=df,
            numeric_cols=numeric_cols,
            categorical_cols=categorical_cols,
            boolean_cols=boolean_cols,
            **kwargs
        )

    if method == 'dbscan':
        return run_dbscan(
            df=df,
            numeric_cols=numeric_cols,
            categorical_cols=categorical_cols,
            boolean_cols=boolean_cols,
            **kwargs
        )

    if method == 'hdbscan':
        return run_hdbscan(
            df=df,
            numeric_cols=numeric_cols,
            categorical_cols=categorical_cols,
            boolean_cols=boolean_cols,
            **kwargs
        )

    if method == 'kprototypes':
        return run_kprototypes(
            df=df,
            numeric_cols=numeric_cols,
            categorical_cols=categorical_cols,
            boolean_cols=boolean_cols,
            **kwargs
        )

    raise ValueError(
        f"Unknown clustering method '{method}'. "
        f"Valid methods: kmeans, agglomerative, dbscan, hdbscan, kprototypes"
    )



def summarize_clusters(
    df: pd.DataFrame,
    cluster_col: str = 'MLCluster'
) -> pd.DataFrame:
    if cluster_col not in df.columns:
        raise ValueError(f"Column '{cluster_col}' not found in DataFrame.")

    summary_rows = []

    for cluster_id, subdf in df.groupby(cluster_col):
        row = {
            'Cluster': cluster_id,
            'Count': len(subdf),
        }

        if 'AW' in subdf.columns:
            row['Mean_AW'] = subdf['AW'].mean()

        if 'Pow' in subdf.columns:
            row['Mean_Pow'] = subdf['Pow'].mean()

        if 'x' in subdf.columns:
            row['Mean_x'] = subdf['x'].mean()

        if 'y' in subdf.columns:
            row['Mean_y'] = subdf['y'].mean()

        if 'AtomName' in subdf.columns:
            row['TopAtomNames'] = ', '.join(
                subdf['AtomName'].value_counts().head(5).index.astype(str).tolist()
            )

        if 'ResidueName' in subdf.columns:
            row['TopResidues'] = ', '.join(
                subdf['ResidueName'].value_counts().head(5).index.astype(str).tolist()
            )

        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)

    if 'Count' in summary_df.columns:
        summary_df = summary_df.sort_values('Count', ascending=False).reset_index(drop=True)

    return summary_df



def compare_manual_vs_ml(
    df: pd.DataFrame,
    manual_col: str = 'CanonicalName',
    ml_col: str = 'MLCluster'
) -> pd.DataFrame:
    ensure_required_columns(df, [manual_col, ml_col])

    return pd.crosstab(df[manual_col], df[ml_col])



def print_cluster_summary(
    summary_df: pd.DataFrame
):
    print("\n=== ML CLUSTER SUMMARY ===")

    if len(summary_df) == 0:
        print("No clusters found.")
        return

    print(summary_df.to_string(index=False))