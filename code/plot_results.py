import matplotlib.pyplot as plt
import pandas as pd


def main():
    filename_ArCI = ".\\results\\ArCI\\ArCI_directionality_results.csv"
    filename_H2O = ".\\results\\H2O\\H2O_directionality_results.csv"

    df_ArCI = pd.read_csv(filename_ArCI)
    df_H2O = pd.read_csv(filename_H2O)

    plots_per_central_group(dfs=[df_ArCI, df_H2O])
    plots_80(df_ArCI, df_H2O)


def plots_80(df1, df2):
    fig, ax1 = plt.subplots(figsize=(10, 5))

    plottitle = "Bin fraction containing 80% data\nAbsolute difference ArCI - H2O"
    plt.title(plottitle)
    plt.xlabel("Resolution")
    ax1.set_ylabel("Fraction of bins")

    contact_groups = df1.contactgroup.unique()

    for contact_group in contact_groups:
        part_df1 = df1[(df1.contactgroup == contact_group)]
        part_df2 = df2[(df2.contactgroup == contact_group)]

        abs_diff = part_df1.bins_80/part_df1.bins - part_df2.bins_80/part_df2.bins

        ax1.plot(part_df1.resolution, abs_diff, label=contact_group)

    ax1.hlines(0, 0.2, 1.5, color="black")
    ax1.legend()
    plt.savefig("bins_80%_data.png")


def plots_per_central_group(dfs):
    for df in dfs:
        contact_groups = df.contactgroup.unique()
        central_groups = df.centralgroup.unique()

        for central_group in central_groups:
            for contact_group in contact_groups:
                plottitle = central_group + " " + contact_group + " stuff"

                part_df = df[(df.centralgroup == central_group) & (df.contactgroup == contact_group)]

                fig, ax1 = plt.subplots(figsize=(10, 5))
                ax2 = ax1.twinx()

                plt.title(plottitle)
                plt.xlabel("Resolution")
                ax1.set_ylabel("Fraction of bins")
                ax2.set_ylabel("Number of bins")

                ax1.plot(part_df.resolution, part_df.bins_in_vdw/part_df.bins, label='Fraction bins in vdw')
                
                ax2.plot(part_df.resolution, part_df.bins_in_vdw*part_df.resolution**3, label='Volume bins in vdw')

                ax1.plot(part_df.resolution, part_df.empty_bins/part_df.bins, label='Fraction empty bins')
                ax1.plot(part_df.resolution, part_df.bins_80/part_df.bins, 
                         label="Fraction bins containing 80% of the data")

                # plt.plot(part_df.resolution, (part_df.bins - part_df.empty_bins)*part_df.resolution**3,
                #          label="%data (1) / volume cluster")

                ax2.plot(part_df.resolution, part_df.datapoints/part_df.bins, label='Nd/Nb', color="purple")
                # ax2.plot(part_df.resolution, part_df.bins, label="Total bins", color="red")

                ax1.legend(loc="upper left")
                ax2.legend(loc="upper right")
                plt.savefig(central_group + "_" + contact_group + ".png")

                print(central_group + " " + contact_group + u'\u2713')


if __name__ == "__main__":
    main()
