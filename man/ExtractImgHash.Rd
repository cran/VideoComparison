\name{ExtractImgHash}
\alias{ExtractImgHash}
\title{
NoSQL JSON oriented database video frame hash extraction function
}
\description{
A R function allowing to ask to a REST/Json server for 
the video frame hashes stored there. 
It will return it a list with the hashes related to the requested image
inside the video for those frames.

Following packages are required: RJSON, RCurl
}
\usage{
ExtractImgHash(pos,father,url="http://localhost:9200/selected_db/selected_db/_search")
}
\arguments{
\item{pos}{String describing the ID of the frame we are interested
into its hashes}
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
\keyword{Image Hash}
\references{
Espinosa-Gutiez, S., Ordieres-Mere, J., Bello-Garcia, A.: Large scale part-to-part video matching by a 
likelihood function using featured based video representation. TAMIDA 2013 - Taller de Mineria de Datos dentro del IV 
congreso Espannol de Informatica. 254-257 (2013). \url{http://bioinspired.dacya.ucm.es/maeb2013/images/ActasCAEPIA_final.pdf}
\cr
Cui, M., Femiani, J., Hu, J., Wonka, P., Razdan, A.: 
Curve matching for open 2D curves. 
Pattern Recogn. Lett. 30, 1-10 (2009)
}
\examples{
%% ~~ simple examples of the most important functions ~~
## Requires specific server architecture
# mm<-ExtractImgHash("000060","C0031D0")
}
