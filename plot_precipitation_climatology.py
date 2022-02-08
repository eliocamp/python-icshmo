import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import argparse


def convert_pr_units(darray):
    """Convert kg m-2 s-1 to mm day-1.

    Args:
      darray (xarray.DataArray): Precipitation data

    """

    darray.data = darray.data * 86400
    darray.attrs['units'] = 'mm/day'

    return darray


def create_plot(clim, model, season, gridlines=False, levels = 1, bins = None):
    """Plot the precipitation climatology.

    Args:
      clim (xarray.DataArray): Precipitation climatology data
      model (str): Name of the climate model
      season (str): Season

    Kwargs:
      gridlines (bool): Select whether to plot gridlines

    """

    fig = plt.figure(figsize=[12,5])
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree(central_longitude=180))

    if bins:
        levels = bins

    clim.sel(season=season).plot.contourf(ax = ax,
                                          levels = levels,
                                          extend='max',
                                          transform = ccrs.PlateCarree(),
                                          cbar_kwargs = {'label': clim.units},
                                          cmap=cmocean.cm.haline_r)
    ax.coastlines()
    if gridlines:
        plt.gca().gridlines()

    title = f'{model} precipitation climatology ({season})'
    plt.title(title)


def main(inargs):
    """Run the program."""

    dset = xr.open_dataset(inargs.pr_file)

    clim = dset['pr'].groupby('time.season').mean('time', keep_attrs=True)
    clim = convert_pr_units(clim)

    create_plot(clim, dset.attrs['source_id'], inargs.season,
                gridlines = inargs.grid, levels = inargs.levels, bins = inargs.bins)
    plt.savefig(inargs.output_file, dpi=200)


if __name__ == '__main__':
    description='Plot the precipitation climatology for a given season.'
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("pr_file", type=str, help="Precipitation data file")
    seasons = ["DJF", "MAM", "JJA", "SON"]
    parser.add_argument("season", type=str, help="Season to plot", choices = seasons)
    parser.add_argument("output_file", type=str, help="Output file name")
    parser.add_argument("-g", "--grid", help = "Would you like gridlines?", action = "store_true")
    parser.add_argument("-l", "--levels", type = float, nargs = "*",
                        default = np.arange(0, 13.5, 1.5),
                        help = "A list of levels to contour.")

    parser.add_argument("-b", "--bins", type = int,
                        help = "Number of bins. Overwrites --levels if not None.",
                        default = None)

    args = parser.parse_args()

    main(args)

