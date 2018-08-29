# em-followup
EM followup of Gravitational Wave sources

# Help

### Module Requirements in Python
Astropy 3.0.3
Requests 2.14.2

### Request_LCO
- This module contains a class for requesting observations with the [Las Cumbres Observatories](https://lco.global) network of telescopes.
- The proposal used in this porject only has access to the 0.4m telescopes, therefore this code is specific to those telescopes.

**Example of use:**

### download_data_LCO
- module to download all LCO images taken with your proposal ID in a time frame
- create txt file called userdata.txt
- write:
- LCOGT archive username, password, datafolder, and the name of the proposals
- datafolder is the name of where you want to save files 
- names of the proposals separated by commas 

- e.g.
- username = 
- password = 
- datafolder = 
- proposals = 

**Example of use:**

```
download(path_general,sdate,edate,proposalID,datafolder)
```

### AAVSO_get_variables
- The Zooniverse project asks users to circle the variations observed in images taken 5-7 days apart.
- This module searches the AAVSO database using given RA and Dec limits and returns a list of coordinates for the known variable stars in 
that field.
- Circle are drawn around these known variables so that only "new" variables will be highlighted in the Zooniverse project.

**Example of use:**
```python
>>>fromRA = 9.1
>>>toRA = 10.8
>>>fromDec = -4.1 	
>>>toDec = 3.2

>>>RA, Dec = get_variables(fromRA, toRA, fromDec, toDec)

>>>print(RA)

9.11692
9.22075
9.32167
9.39625
9.45395

>>>print(Dec)

-0.52356
3.11739
1.96000
0.29417
-0.05333
```

### imcreate
- this module processes images to make them compatible and produces difference images, saving them as .png files.
- create a general folder that contains the following empty folders:
    
- 'images'
- 'LCO_images'
- 'fits_images'
- 'subtraction'
- 'relative'
- 'Zooniverse_upload'
- 'rejects'
- 'finished_images'

- In 'images', create a folder per target to store target images
- 2 comparison images of single target per target folder
- images must have same name keyword as folder
- e.g. target folder called 'keyword', images within called 'keyword_1', 'keyword_2'
- name target folders using actual target name 


**Example of use:**
```
path_general='/Users/lewisprole/Documents/University/year3/summer_project'
remove_space(path_general,'images')
remove_space(path_general,'LCO_images/raw')
combine_fold(path_general)
reject_dir(path_general)
fz_remove(path_general)
image_rename(path_general)
CPR(path_general)
bright_diff(path_general)
sub(path_general)
out_save(path_general,'yes','yes')
```
### run_code

- script to run the various python scripts to request, download and process images
- all you need to do is enter your details and run the script 

```
#requests details
username= "username"
password = "password"
RA = 202.4708
Dec = 47.1953
PROPOSAL_ID = "proposalID"
magnitude = 9
expt=30
period =10000 #how often to check the if the images have downloaded
path_general='path/to/general_folder'

#download details
sdate='2018-06-19'
edate='2018-08-19'
proposalID="LCOEPO2018A-004"
datafolder="path/to/LCO_images"
```