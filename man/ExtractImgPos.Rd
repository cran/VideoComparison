\name{ExtractImgPos}
\alias{ExtractImgPos}
\title{
NoSQL JSON oriented database video frame extraction function
}
\description{
A R function allowing to ask to a REST/Json server for 
the video frame descriptions stored there. 
It will return it as a vector of strings making clear the positions
inside the video for those frames.

Following packages are required: RJSON, RCurl
}
\usage{
ExtractImgPos(father,url="http://localhost:9200/selected_db/selected_db/_search")
}
\arguments{
\item{father}{String describing the ID key allowing to access the data
into the database.}
\item{url}{String describing the URL which will answer queries in JSON 
oriented mode.}
}
\author{
Joaquin Ordieres, <j.ordieres@upm.es>,
Antonio Bello, <abello@uniovi.es>.
}
\keyword{JSON}
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
## Requires specific server architecture
# mm<-ExtractImgPos("C0031D0")
}
