import seaborn as sns
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.figure as mfigure
import numpy as np
import numpy as np
import pandas as pd
from typing import Optional
from datetime import datetime
import polars as pl
from fastapi import HTTPException
from mhw_utils import process_mhw_data
from io import BytesIO
# from tempfile import NamedTemporaryFile


def convert_to_day(date_str: str) -> str:
    try:
        datetime_obj = datetime.strptime(date_str, '%Y-%m-%d')
        return date_str
    except ValueError:
        if len(date_str) == 4:
            format_str = '%Y'
        elif len(date_str) == 6:
            format_str = '%Y%m'
        elif len(date_str) == 8:
            format_str = '%Y%m%d'
        else:
            raise ValueError("Invalid date format")

        # Convert the date string to a date object
        try:
            date_obj = datetime.strptime(date_str, format_str)
            # Return the date in the 'YYYY-MM-DD' format
            return date_obj.strftime('%Y-%m-%d')
        except ValueError as e:
            raise ValueError("Invalid date: " + str(e))


def period2date(period: str, init: str = '1982-01-01'):
    start, end = period.split('-')
    try:
        if not start:
            start = init
        else:
            start = convert_to_day(start)

        if not end:
            end = str((datetime.today() - pd.DateOffset(months=1)
                       ).replace(day=1).date())
        else:
            end = convert_to_day(end)

        return start, end

    except ValueError as e:
        raise HTTPException(
            status_code=400, detail=str(e))


def genImageFile(fig):
    # Save plot to a bytes buffer
    buf = BytesIO()
    fig.savefig(buf, format='png')
    # Yes, it's required to reset the cursor to the beginning of the buffer.
    buf.seek(0)
    # Write the image data to a temporary file
    # with NamedTemporaryFile(delete=False, suffix=".png") as temp_file:
    #    temp_file.write(buf.read())
    #    temp_file_name = temp_file.name
    # Return the image as a file response
    # return FileResponse(temp_file_name, media_type="image/png")
    return buf


async def month_climatology(lon0: float, lat0: float, lon1: Optional[float] = None, lat1: Optional[float] = None, period: Optional[str] = "all", sstvar: Optional[str] = None, palette: Optional[str] = None):
    # Get data
    enddate = (datetime.today() - pd.DateOffset(months=1)
               ).replace(day=1).date()
    if not sstvar or sstvar.strip() not in ['sst', 'sst_anomaly']:
        sstvar = 'sst'

    _, data = await process_mhw_data(lon0=lon0, lat0=lat0, lon1=lon1, lat1=lat1,
                                     start="1982-01-01", end=str(enddate), append=sstvar, mode="month_mean")

    # data = data.with_columns(pl.col('date').str.strptime(pl.Date, "%Y-%m-%dT%H:%M:%S", strict=False))
    data = data.with_columns([
        pl.col('date').dt.year().alias("year"),
        pl.col('date').dt.month().alias("month")
    ])

    if not period or period == "all":
        period = [f"{data['date'].min().year}-{data['date'].max().year}"]
    else:
        # [var.strip() for var in period.split(',')] #strip in loop
        period = period.split(',')

    if not palette:
        palette = sns.cubehelix_palette(len(period), gamma=.5)
    else:
        palette = [var.strip() for var in palette.split(',')]
        if len(palette) < len(period):
            palette = palette + \
                sns.cubehelix_palette(len(period) - len(palette), gamma=.5)
        elif len(palette) > len(period):
            palette = palette[0:len(period)]

    fig, ax = plt.subplots(figsize=(12, 8))

    for i, per in enumerate(period):
        per = per.strip()
        chkyear = True
        try:
            int(per)
        except ValueError:
            chkyear = False

        if chkyear:
            # Current year, plot with red
            yr = int(per[0:4])
            subset = data.filter(data['year'] == yr)
            color = "red" if yr == enddate.year else palette[i]
            label = str(per)
        else:
            start_year, end_year = map(lambda x: int(x[0:4]), per.split("-"))
            subset = data.filter(
                (data['year'] >= start_year) & (data['year'] <= end_year))
            color = palette[i]
            label = f"{start_year}-{end_year}"

        if not subset.is_empty():
            monthly_avg = subset.groupby('month').agg(
                pl.col(sstvar).mean()).sort('month')
            ax.plot(monthly_avg['month'],
                    monthly_avg[sstvar], color=color, label=label)

    if sstvar == 'sst':
        ylab = 'SST'
    else:
        ylab = 'SST Anomaly'

    ax.legend()
    ax.set_xlabel('Month')
    ax.set_ylabel('Average ' + ylab)
    ax.set_title(
        f"Average Monthly {ylab} at Lon: {str(lon0) if not lon1 else str(lon0) + ' - ' + str(lon1)}, Lat: {str(lat0) if not lat1 else str(lat0) + ' - ' + str(lat1)}")
    ax.set_xticks(np.arange(1, 13))
    ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr',
                        'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    return genImageFile(fig)


async def region_climatology(bbox: str, start: str, end: str, sstvar: Optional[str] = None, palette: Optional[str] = None):
    bbox_list = bbox.strip('()').split('),(')
    if not sstvar or sstvar.strip() not in ['sst', 'sst_anomaly']:
        sstvar = 'sst'

    if not palette:
        palette = sns.cubehelix_palette(len(bbox_list), gamma=.5)
    else:
        palette = [var.strip() for var in palette.split(',')]
        if len(palette) < len(bbox_list):
            palette = palette + \
                sns.cubehelix_palette(len(bbox_list) - len(palette), gamma=.5)
        elif len(palette) > len(bbox_list):
            palette = palette[0:len(bbox_list)]

    # fig, ax = plt.subplots(figsize=(12, 8))
    fig = mfigure.Figure(figsize=(12, 8))
    ax = fig.subplots()

    for i, bbox_str in enumerate(bbox_list):
        bbox_coords = [float(coord) for coord in bbox_str.split(',')]
        if len(bbox_coords) == 2:
            lon0, lat0, lon1, lat1 = bbox_coords[0], bbox_coords[1], None, None
        elif len(bbox_coords) == 4:
            lon0, lat0, lon1, lat1 = bbox_coords
        else:
            raise HTTPException(
                status_code=400, detail="Invalid bounding box format. Must contain 2 or 4 coordinates.")

        # print(lon0, lon1, lat0, lat1, start, end)
        _, df = await process_mhw_data(lon0=lon0, lon1=lon1, lat0=lat0, lat1=lat1,
                                       start=start, end=end, append=sstvar, mode='area_mean_'+sstvar)
        df = df.sort('date')  # Ensure data is in chronological order

        # Convert date to year-month format for plotting
        # df = df.with_columns(pl.col('date').dt.strftime('%Y-%m').alias('year_month'))
        ax.plot(df['date'], df[sstvar].to_list(),
                color=palette[i], label=f'BBOX{i+1}: ({bbox_str})')

    if sstvar == 'sst':
        ylab = 'SST'
    else:
        ylab = 'SST Anomaly'

    # Set the date formatter for the x-axis to only display the year
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    # Set the major ticks to be every 5 years
    # ax.xaxis.set_major_locator(mdates.YearLocator(5))
    major_ticks = mdates.YearLocator()
    ax.xaxis.set_major_locator(major_ticks)

    # Set the minor ticks to be each month
    ax.xaxis.set_minor_locator(mdates.MonthLocator())

    # Increase the length of the major ticks
    ax.tick_params(which='major', length=10)

    # Only display labels for every 5th major tick
    for i, label in enumerate(ax.xaxis.get_majorticklabels()):
        if int(label.get_text()) % 5 != 0:
            label.set_visible(False)

    ax.tick_params(which='major', length=10)
    ax.set_xlabel("Year")
    ax.set_ylabel(ylab)
    ax.set_title(f"{ylab} Climatology from {start} to {end}")
    fig.autofmt_xdate()
    # plt.show()
    return genImageFile(fig)
