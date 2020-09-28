# gmslCMIP6
This repository can be used for pulling CMIP6 data from ESGF, preprocessing raw (monthly mean, but probably also other frequencies) 'zostoga' and 'tas' data, including time-merging, dedrifting and calculating area-weighted means, and analyzing the preprocessed data.

Instructions for running the WGET script can be found here: https://esgf.github.io/esgf-user-support/user_guide.html#download-data-from-esgf-using-wget, for example using:

```
bash wget-20200715133532_zostoga.sh
```

will pull zostoga files from ESGF and store them per variable per model. Note that a 5 second delay has been built in between each download, to avoid overasking ESGF nodes. The script may require a few iterations before the download is complete, depending on server response.

Required Python packages: xarray, numpy, fnmatch, cdo
