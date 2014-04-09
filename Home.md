#Installation instructions for metaphlan in a galaxy environment.


Clone the metaphlan galaxy repository 
```
https://bitbucket.org/george_weingart/metaphlan_galaxy 
```
into the ~/galaxy-dist/tools directory

Rename the resulting metaphlan_galaxy directory to "metaphlan"

Update member tool_conf.xml  in the galaxy directory adding the following: 

```
    <section name="MetaPhlAn" id="metaphlan">
		<tool file="metaphlan/metaphlan.xml"/>
    </section>
```


Recycle galaxy




Note 

****

bowtie2 is a strict requirement - if missing, please install it:

sudo apt-get install bowtie2


