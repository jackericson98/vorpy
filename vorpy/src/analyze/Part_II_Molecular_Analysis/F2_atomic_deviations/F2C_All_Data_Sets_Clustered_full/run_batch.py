import os

from .config import DATA_ROOT, MOLECULE_GROUPS, OUTPUT_ROOT, SETTINGS_COMBINATIONS
from .pipeline import make_settings_name, main, resolve_group_folders


def run_batch_cluster_plots():
    os.makedirs(OUTPUT_ROOT, exist_ok=True)

    for molecule_label, rel_paths in MOLECULE_GROUPS.items():
        folders = resolve_group_folders(DATA_ROOT, rel_paths)

        if len(folders) == 0:
            print(f"Skipping {molecule_label}: no folders found.")
            continue

        print(f"\n=== RUNNING GROUP: {molecule_label} ===")
        print(f"Folders: {folders}")

        for cfg in SETTINGS_COMBINATIONS:
            mode = cfg['mode']
            method = cfg['ml_method']

            out_subdir = os.path.join(OUTPUT_ROOT, molecule_label, method)
            os.makedirs(out_subdir, exist_ok=True)

            settings_name = make_settings_name(
                mode=mode,
                ml_method=method,
                n_clusters=cfg['n_clusters'],
                eps=cfg.get('eps'),
                min_samples=cfg.get('min_samples'),
                min_cluster_size=cfg.get('min_cluster_size'),
                numeric_cols=cfg['numerical_cols'],
                categorical_cols=cfg['categorical_cols'],
                use_sol_binary=cfg['use_sol_binary'],
                point_alpha=0.2,
                ellipse_n_std=1.2,
            )

            print(f"\nRunning: {molecule_label} | {settings_name}")

            try:
                main(
                    folders=folders,
                    atom_name_field='Name',
                    molecule_class='small_molecule' if molecule_label == 'small molecule' else (
                        'dna' if molecule_label == 'dna' else (
                            'rna' if molecule_label == 'rna' else 'protein'
                        )
                    ),
                    volume_range=(3, 22),
                    save_csv=False,
                    save_plot=True,
                    show_points=True,
                    show_numbers=True,
                    annotate_direct_groups=False,
                    plot_min_count=50,
                    max_spread=3.5,
                    ellipse_min_count=50,
                    ellipse_max_spread=None,
                    ellipse_n_std=1.2,
                    point_alpha=0.2,
                    use_ml_clustering=cfg['use_ml_clustering'],
                    ml_method=method if method != 'manual' else 'kmeans',
                    use_sol_binary=cfg['use_sol_binary'],
                    n_clusters=cfg['n_clusters'] if cfg['n_clusters'] is not None else 12,
                    min_samples=cfg.get('min_samples', 10) if cfg.get('min_samples') is not None else 10,
                    eps=cfg.get('eps', 1.0) if cfg.get('eps') is not None else 1.0,
                    min_cluster_size=cfg.get('min_cluster_size', 25) if cfg.get('min_cluster_size') is not None else 25,
                    numerical_cols=cfg['numerical_cols'],
                    categorical_cols=cfg['categorical_cols'],
                    boolean_cols=cfg['boolean_cols'],
                    output_base=out_subdir,
                    output_name=settings_name,
                )
            except Exception as e:
                print(f"\nFAILED: {molecule_label} | {settings_name}")
                print(f"Reason: {e}")

                fail_log = os.path.join(out_subdir, "_failed_runs.txt")
                with open(fail_log, "a", encoding="utf-8") as f:
                    f.write(f"{molecule_label} | {settings_name}\n")
                    f.write(f"{type(e).__name__}: {e}\n\n")

                continue


if __name__ == "__main__":
    run_batch_cluster_plots()
