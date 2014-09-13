import sys

import userinterface as ui
from acquisitionimage import AcquisitionImage
from selection import SelectionCursorListener

def main(argv=[__name__]):
    acqimage = AcquisitionImage(argv[1])

    ui.display(acqimage.get_science_data(), zscale=True)

    # the listener will collect the actual data we want
    verbose = False
    listener = SelectionCursorListener(acqimage, verbose)
    listener.start()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
