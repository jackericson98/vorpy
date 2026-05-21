from .pipeline import main


if __name__ == "__main__":
    main(
        folders=None,
        atom_name_field='Name',
        molecule_class='dna',
        volume_range=(3, 22),
        save_csv=True,
        save_plot=True,
        show_points=True,
        show_numbers=True,
        annotate_direct_groups=False,
        plot_min_count=50,
        max_spread=3.5,
        ellipse_min_count=50,
        ellipse_max_spread=None,
        ellipse_n_std=2,
        point_alpha=0.2,
        use_ml_clustering=True,
        ml_method='kmeans',
        use_sol_binary=False,
        n_clusters=12,
        min_samples=10,
        eps=1.0,
        min_cluster_size=25,
        numerical_cols=['AW', 'Pow', 'DeltaV'],
        categorical_cols=[],
        boolean_cols=[],
    )
