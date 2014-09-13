
Installation
============

This is a preliminary release of the NICI Reduction package.

The NICI module comes as a tar file: nici.tar.gz


1. Grab the tar file into your local directory.

2. Untar 

 ::

  tar xf nici.tar.gz

3. Setting up the scripts:

 - Define the unix environment variable *nicipath* that points to the directory
   where you untar  nici and add this path to your PYTHONPATH. Please put   
   these in your .cshrc if you can.

 ::       

     (edit .cshrc)

       setenv nicipath <nicipath>     
       setenv PYTHONPATH {$PYTHONPATH}:{$nicipath}

     Example: You untared in /home/user/gem_python
       setenv nicipath /home/user/gem_python
       setenv PYTHONPATH {$PYTHONPATH}:{$nicipath}
      
 - Now update your environment variables:

 ::       

       source .cshrc

 - Define the reduction scripts by sourcing *nici.csh*:

 ::       

       cd $nicipath
       source nici.csh

4. This will define aliases to the nici scripts:

 - alias ncqlook    $nicipath/nici/ncqlook.py
 - alias ncprepare  $nicipath/nici/ncprepare.py
 - alias ncmkflats  $nicipath/nici/ncmkflats.py
 - alias ncscience  $nicipath/nici/ncscience.py

5. Before running the nici scripts.

   Make sure that your can run python specially if you
   can access numpy and scipy packages.

6. For help please see the URL:

 ::       

     file://$nicipath/nici/doc/build/html/index.html  

7. Test:  ncqlook -h

   Should give you a menu of the scrip parameters.

8. PROBLEMS and support.

   -- Make sure you have an updated version of Python with scipy and numpy running.

   -- If test7 fails with a Python problem; e.g. library not found, most probably 
   is due to an incomplete or incompatible installation of the above 2 modules.

   -- For Python support please see:
      http://www.stsci.edu/resources/software_hardware/pyraf/stsci_python

   -- For the nici scripts support please contact: nzarate@gemini.edu

