Installation instructions for metaphlan in a galaxy environment.
These instructions require the Mercurial versioning system, galaxy, and an internet connection.

1.  Clone the metaphlan galaxy repository (https://bitbucket.org/george_weingart/metaphlan_galaxy) into the ~/galaxy-dist/tools directory 
2.  Rename the resulting metaphlan_galaxy directory to "metaphlan"

3. Update member tool_conf.xml  in the galaxy directory adding the following: 

    <section name="MetaPhlAn" id="metaphlan">
		<tool file="metaphlan/metaphlan.xml"/>
	</section>

4. Recycle galaxy

Note #1
******* 
The bowtie2 databases must be posted in the directory ~/galaxy/tools/metaphlan
~/galaxy-dist/tools/metaphlan/bowtie2db$ ls
mpa.1.bt2  mpa.2.bt2  mpa.3.bt2  mpa.4.bt2  mpa.rev.1.bt2  mpa.rev.2.bt2


Note #2
*******
bowtie2 is a strict requirement - if missing, please install it:
sudo apt-get install bowtie2
We need to resolve the issue of how to post the databases in bitbucket
