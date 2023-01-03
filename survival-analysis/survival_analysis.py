from lifelines import CoxPHFitter
from lifelines.datasets import load_rossi


if __name__ == '__main__':
    data = load_rossi()
    cph = CoxPHFitter()
    cph.fit(data, duration_col='week', event_col='arrest')
    cph.print_summary()