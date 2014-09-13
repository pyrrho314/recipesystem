  Nelson,

If we need such an option for testing or debugging, then I'm happy with it being there, but I don't think we should limit the number of detections for general use. The catalogs it generates will get stored with the file and potentially used later on, so they should be complete rather than arbitrarily cut off at some source count.

Of course, it could be argued that the algorithm itself could contain some kind of sanity check - if it starts detecting more sources than there are pixels, then something has clearly gone wrong for example... :-) If we need to set some sanity cut-off, then I'm OK with that, but we should be careful to make sure that it's significantly higher than the maximum number of real sources that could be detected.

Of course, source / noise differentiation can be difficult problem, especially as you point out in crowded fields. I've no real experience with daofind et al, I've always used sextractor. with that experience in mind, I would be optimistic that we could find a good solution. I would imagine it might be necessary to generate statistics from the image to use in determining parameters passed to the source detection task, and also that determining the parameters for the source detection algoritm might be different for different wavelength regimes and instruments...

cheers,
Paul.



On 12/03/2010 09:27 AM, Nelson Zarate wrote:
> Hi Paul;
>
> Do you want to put a selection criteria on a maximum number of sources that
> are found in a given image?
> I am thinking on a crowded field where speckle just above the background
> could
> be found as objects. I am somewhat familiar with daofind but we can use
> some
> other like sextractor.
>
> The other question is about flux. Daofind gives magnitude. What would be
> the
> transformation?
>
> Thanks
> Nelson
Nelson,

Yes, putting fwhm in the object table alongside flux is a good plan.

With regard to calculating a "Seeing" number...:

The main thing to be aware of is that the current code the DAS have I think has the cataloging step rigged (presumably by setting parameters appropriately) to only detect star like objects and to reject galaxies and other non point sources.

Now of course the proper(tm) way to do this is to detect all the sources, and then to do something more sophisticated than an arithmetic mean when doing the average. Essentially, the histogram of FWHMs will normally look like a very skewed gausian, with a long tail to high fwhm. - the main peak is the point sources, then the tail is the galaxies and other extended sources.

So, I'm not sure off hand what the best solution is. Finding the modal value might be a good approach, or fitting a gaussian to the histogram, or more pragmatically, a (multiply) clipped mean might be good, but I've not really looked into this recently - it may be that there's a good accepted practice that works well or we may have to figure it out for ourselves...

Either way, I'd hold off on putting in seeing as a simple arithmetic mean for now. Probably the best way forward would be to get it working with the fwhm in the object table, then we should generate the object tables for say a dozen or so different gmos frames, calculate the above statistics on their fwhm values and compare them to what the DASs or scientists say the seeing on those frames is (picking stellar like objects by eye), and choose the statistic that matches best. :-)

Extragalactic fields, especially with galaxy clusters are going to be the most difficult to get right...

to answer the specific question, no, we don't need the reference catalog for this.

Cheers,
 Paul.



Hi Nelson,

Here's an outline of how I think it should work:

Cheers,
 Paul.



-------- Original Message --------
Subject: Quick bullet point thoughts re flux calibration
Date: Fri, 19 Nov 2010 12:12:01 -1000
From: Paul Hirst <phirst@gemini.edu>
To: <klabrie@gemini.edu>

Hi Kathleen,

Here's some quick thoughts as to how I think we should be laying this
out:
I do realize that this will result in a lot of extra fits extensions
with catalog tables in them. However, I don't think that's a problem...

At the recipe level, something like:

# All the usual single frame processing
# including flat field etc
# Should have an approximate (from the telescope WCS) too at this point

detect_sources
# This would try to detect all the sources in the frame and would
generate a source catalog that should be inserted into the output file
as a fits table of some kind, with one fits table catalog extension for
each image extension. The catalog should contain at least the following:
# source ID number (simple serial number unique within this table)
# Centroid X pixel co-ordinate
# Centroid Y pixel co-ordinate
# Source flux (in ADUs / whatever units the data are in at this point)
# error on flux measurement
# If the cataloger measures flux using object and sky appertures, then
these values should probably be listed too
# RA (assume WCS in image is correct)
# Dec (ditto)
# Any kind of detection quality flags that the cataloger produces

measure_iq
# In the long term, the measure IQ primitive should simply use this
catalog as input. In the even longer term, this might move later in the
recipe to take advantage of the reference catalog data.

add_reference_catalogs
# Initially, the only reference catalog should be the gmos flux
standards list. Later we will add SDSS, and possibly others.
# For each image extension, this primitive should simply query each/the
reference catalog for all sources that fall within the footprint of the
image (plus some margin to allow for a dodgy WCS), and simply add a fits
table extension to the data file containing that catalog data.

correlate_with_reference_catalogs
# This would correlate the source detection table with the reference
catalog table. I'm open to suggestions as to the exact format of the
output, but to keep things fairly simple I would suggest that for each
reference catalog, it adds a column to the source detection table, into
which it writes the reference catalog ID that it thinks corresponds to
that detected source.

calculate_zeropoint
# For each detected source with a reference catalog association,
calculate the zeropoint. Calculate mean and sigma over all detected
sources with associations. These should be added as header keywords to
the image extensions they are derived from, and of course logged etc.

refine_wcs
# Let's not worry about this for now, but at this point one would use
the association between detections and reference sources to refine the
WCS in the image.



Hope this helps,
  Paul.


On 11/29/2010 09:49 AM, Nelson Zarate wrote:
> Hi Paul;
>
> I assumed you received the attach mail from Kathleen.
>
> I had the idea of merging the 2 versions that that GN and GS have
> written for the ZP correction but it seems that you have something
> different about
> this calculation.
>
> Do you have requirements for this?
>
> Thanks
> Nelson
>
>

