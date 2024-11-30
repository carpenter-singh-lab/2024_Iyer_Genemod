import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

hyperparams = ["nhead", "nhid", "nlayers", "dropout", "learning_rate"]
SELECTED_RUN = "795f9503ff2ddc291ac9fb64fc5a36ed"


def plot_loss_grid(df, metric="loss"):
    """Create a grid plot of all training runs

    Args:
        df: DataFrame containing the training data
        metric: String indicating which metric to plot (default: 'loss')
    """

    # Get unique folders
    folders = df["folder"].unique()
    n_runs = len(folders)

    # Calculate grid dimensions (trying to make it roughly square)
    n_cols = int(np.ceil(np.sqrt(n_runs)))
    n_rows = int(np.ceil(n_runs / n_cols))

    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 20), sharey=True, sharex=True)

    # Flatten axes for easier iteration if multiple rows
    axes_flat = axes.flatten() if n_rows > 1 else axes

    # Define short names for hyperparams
    param_abbrev = {
        "nhead": "nh",
        "nhid": "hid",
        "nlayers": "nl",
        "dropout": "dr",
        "learning_rate": "lr",
    }

    # Plot each run
    for idx, folder in enumerate(folders):
        if idx < len(axes_flat):
            ax = axes_flat[idx]
            run_data = df[df["folder"] == folder]

            # Get hyperparams for this run
            params = {k: run_data[k].iloc[0] for k in hyperparams}
            param_str = ",".join(
                [f"{param_abbrev[k]}:{params[k]}" for k in hyperparams]
            )

            # Plot specified metric
            assert (
                run_data["val_loss"] == run_data["val_loss"].iloc[0]
            ).all(), "Final validation loss is not the same across the whole run"
            val_score_selected = run_data["val_loss"].iloc[0]

            epoch_selected = run_data[run_data["score"] == val_score_selected][
                "epoch"
            ].iloc[0]

            metric_data = run_data[
                (run_data["metric"] == metric) & (run_data["epoch"] <= epoch_selected)
            ]

            for subset, color in zip(["train", "valid"], ["blue", "red"]):
                subset_data = metric_data[metric_data["subset"] == subset]
                ax.plot(
                    subset_data["epoch"],
                    subset_data["score"],
                    color=color,
                    alpha=0.8,
                    label=subset,
                )

            # Update title to include both run ID and params
            if folder == SELECTED_RUN:
                ax.set_facecolor("yellow")

            ax.set_title(f"Run {folder[:8]}\n{param_str}", fontsize=6)
            ax.tick_params(labelsize=6)
            ax.grid(True, alpha=0.3)

            if idx == 0:
                ax.legend(fontsize=6)

    # Remove unused subplots
    for j in range(len(folders), len(axes_flat)):
        fig.delaxes(axes_flat[j])

    # Add common labels
    fig.text(0.5, 0.0, "Epoch", ha="center", fontsize=12)
    fig.text(0.0, 0.5, metric, va="center", rotation="vertical", fontsize=12)

    plt.tight_layout()
    import re

    metric_ = re.sub(r"\W+", "_", metric)
    plt.savefig(f"output/grid_{metric_}.png", dpi=300, bbox_inches="tight")
    plt.close()


def main():
    # Read data
    from pathlib import Path

    df_learning = pd.read_csv(
        Path(__file__).parent / "output/leave_compound_out_learning_log.csv.gz"
    )
    df_learning["metric"] = df_learning["metric"].str.replace("@top{k}", "_top_k")
    df_configs = pd.read_csv(
        Path(__file__).parent / "output/leave_compound_out_combined.csv"
    )

    df_configs = df_configs.rename(columns={"hash_id": "folder"})
    df = df_learning.merge(df_configs, on="folder")

    # specify the sort order of the rows, based on a list of columns
    df = df.sort_values(by=hyperparams)

    # ensure that the number of rows matches
    assert len(df) == len(df_learning)

    # Generate grids
    plot_loss_grid(df, metric="loss")
    plot_loss_grid(df, metric="auprc")
    plot_loss_grid(df, metric="p_top_k")

    # Analyze overfitting
    # results, summary = analyze_overfitting(df)

    # print("\nFinal Gap Analysis Summary:")
    # print(tabulate(summary, headers=summary.columns, tablefmt="pipe", floatfmt=".4f"))

    # # Save results
    # results.to_csv("output/learning_gap_analysis.csv", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
