\name{VideoComparison}
\alias{VideoComparison}
\title{
Calculate the matching segment between two
video motion vectors according the the minimum length required.
}
\description{
For two video segments represented by their own video motion vector
it is interesting to look for their matching in order to verify
potential coherence between them happens.

The video motion curves can be derived from the opensource 
tool called shotdetect \url{http://johmathe.name/shotdetect.html}.

By using an improved algorithm based on the Cui et al paper
it is implemented the identification of which segment 
bigger than the requested length, is common to both videos, 
in terms of motion structure.
}
\usage{
VideoComparison(mm,m2,stp,nprocesses=0)
}
\arguments{
\item{mm}{Motion vector for video 1. It is supposed to be contained in the second one so, it will be shorter than the second vector}
\item{m2}{Motion vector for video 2. It is supposed to be the container}
\item{stp}{Unit of motion curvature to be considered as relevant 
  for our analysis.}
\item{nprocesses}{Number of processes that should be spawned. If 0, it's left to the function's discretion. Requires package parallel to spawn more than 1 process.}
}
\details{
This function maximizes the score (Score of matching = length of matching * degree of matching), but returns it separated by its components.

If the package parallel is available, this function will run across several processes.
}
\value{
The output provided is a list with following components:
\item{$sc}{Degree of matching}
\item{$pos1}{Frame where matching starts in video1}
\item{$pos2}{Frame where matching starts in video2}
\item{$lngth}{Length of matching segment found}
}
\author{
Silvia Esinosa, <silviaespinosag@gmail.com>,
Jose Maria Perez, <josem.perez.ramos@gmail.com>,
Joaquin Ordieres, <j.ordieres@upm.es>,
Antonio Bello, <abello@uniovi.es>.
}
\keyword{Video Motion}
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

#load example json from data folder
fileName1 = strcat(path.package(package="VideoComparison"),"/extdata/videomotion1.gz")
json1 = readLines(gzfile(fileName1))
fileName2 = strcat(path.package(package="VideoComparison"),"/extdata/videomotion2.gz")
json2 = readLines(gzfile(fileName2))

#Extract data from json
hh<- as.numeric(unlist(lapply(fromJSON(json1)$hits$hits[[1]]$`_source`$video_hash$frames,head,1)))
h2<- as.numeric(unlist(lapply(fromJSON(json2)$hits$hits[[1]]$`_source`$video_hash$frames,head,1)))

## Requires specific server architecture
# hh<-ExtractMotion("C0031D0")
# h2<-ExtractMotion("C0035D0")

mm<-VideoComparison(hh,h2,12,nprocesses=1);

}
