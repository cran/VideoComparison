\name{Compare2Videos}
\alias{Compare2Videos}
\title{
Compare two videos starting from their names.
}
\description{
Compare two videos from their names by retrieving their data 
from the database.
A likelihood estimation will be returned, as well as 
basic parameters for motion matching (correlation, pos on each
video and length of the matching).
}
\usage{
Compare2Videos(nv1,nv2,stp=10,fsc=0)
}
\arguments{
\item{nv1}{String with the identification of the first name.
}
\item{nv2}{String with the identification of the second name.
}
\item{stp}{Unit of motion curvature to be considered.}
\item{fsc}{Way to determine the correlation between videos.
When 0 it will adopted the two times of score derived from
the motion video analysis. Otherwise the specified factor
will be used.
}
}
\author{
Joaquin Ordieres, <j.ordieres@upm.es>,
Antonio Bello, <abello@uniovi.es>.
}
\keyword{NoSQL}
\keyword{Video Likelihood}
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
## Requires specific server architecture
# lv<-VideoSearch()
# Compare2Videos(lv[1],lv[2])
}
