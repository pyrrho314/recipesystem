from testutil import assert_tolerance, get_data_file_name
from astrodata import AstroData

from selection import get_selection_peak

def test_finding_the_center_of_titan():
    ad = AstroData(get_data_file_name("N20060131S0011.fits"))
    selection = get_selection_peak(ad.data, (391.0, 539.0), float(ad.pixel_scale()))

    predicted_center = selection.get_center()

    # the object is titan, doesn't have a very clearly defined center
    assert_tolerance(predicted_center,
                     (384.87060478881705, 542.73266988305159),
                     tolerance=0.07)

def test_finding_an_easier_center():
    ad = AstroData(get_data_file_name("N20060131S0012.fits"))
    selection = get_selection_peak(ad.data, (489.07, 478.99), float(ad.pixel_scale()))

    predicted_center = selection.get_center()

    assert_tolerance(predicted_center,
                     (489.11025833697016, 478.68088198208636),
                     tolerance=0.02)


def test_out_of_bound_aperture():
    ad = AstroData(get_data_file_name("N20060131S0012.fits"))
    selection = get_selection_peak(ad.data, (10.0, 10.0), float(ad.pixel_scale()))

    predicted_center = selection.get_center()

    # should be rather non-sense that is returned
    assert_tolerance(predicted_center,
                     (50.0, 50.0),
                     tolerance=50.0)
