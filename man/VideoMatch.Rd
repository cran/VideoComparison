\name{VideoMatch}
\alias{VideoMatch}
\title{
Starting from two list of image hashes and the score in motion matching 
factor.
}
\description{
VideoMatch looks for the likelihood of coherence between the sequence
of pertinent image hashes obtained after preparing the comparison between
video Motion vector.

Likelihood estimation relates the coherence (distance: norlized average 
of four normalized distances) between image sequence from the two videos
in the region where their motion curve matches. It can be affected by 
the factor of coherence of this matching.
}
\usage{
VideoMatch(lh1,lh2,sc)
}
\arguments{
\item{lh1}{List with the ordered sequence on hashes 
from images belonging to the common motion chunck.
It refers to the video1.
}
\item{lh2}{List with the ordered sequence on hashes 
from images belonging to the common motion chunck.
It refers to the video2.
}
\item{sc}{This is the motion scoring factor.
It can be derived from the VideoComparison function
(first argument of the answer).
}
}
\author{
Joaquin Ordieres, <j.ordieres@upm.es>,
Antonio Bello, <abello@uniovi.es>.
}
\keyword{NoSQL}
\keyword{Video Motion Curve}
\references{
Espinosa-Gutiez, S., Ordieres-Mere, J., Bello-Garcia, A.: Large scale part-to-part video matching by a 
likelihood function using featured based video representation. TAMIDA 2013 - Taller de Mineria de Datos dentro del IV 
congreso Espannol de Informatica. 254-257 (2013). \url{http://bioinspired.dacya.ucm.es/maeb2013/images/ActasCAEPIA_final.pdf}
\cr
Cui, M., Femiani, J., Hu, J., Wonka, P., Razdan, A.: Curve matching for 
open 2D curves. Pattern Recogn. Lett. 30, 1-10 (2009)
}
\examples{
%% ~~ simple examples of the most important functions ~~
#Load json data from txt
fileName1 = strcat(path.package(package="VideoComparison"),"/data/imagehash1.txt.gz")
json1 = readLines(gzfile(fileName1))
fileName2 = strcat(path.package(package="VideoComparison"),"/data/imagehash2.txt.gz")
json2 = readLines(gzfile(fileName2))

# Extract data from json
out1<-fromJSON(json1)$hits$hits[[1]]$`_source`$Hash;
dct1<-as.character(out1[1])
hstrada1<-as.numeric(unlist(strsplit(out1[2],",")))
mw1<-as.numeric(unlist(strsplit(out1[3],",")))
rd1<-as.numeric(unlist(strsplit(out1[4],",")))
img1<-list(dct=dct1,hstrada=hstrada1,mw=mw1,rd=rd1)

out2<-fromJSON(json2)$hits$hits[[1]]$`_source`$Hash;
dct2<-as.character(out2[1])
hstrada2 <-as.numeric(unlist(strsplit(out2[2],",")))
mw2<-as.numeric(unlist(strsplit(out2[3],",")))
rd2<-as.numeric(unlist(strsplit(out2[4],",")))
img2<-list(dct=dct2,hstrada=hstrada2,mw=mw2,rd=rd2)

## Requires specific server architecture
# img1<-ExtractImgHash("000060","C0031D0")
# img2<-ExtractImgHash("000009","C0035D0")

VideoMatch(img1,img2)
}
