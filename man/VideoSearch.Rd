\name{VideoSearch}
\alias{VideoSearch}
\title{
Obtain the ID for main video entries into the database.
}
\description{
Video  hashes have been stored in differnt chuncks. Main chunk
includes the container data. Then, a second chunk named _video.json
stores the motion vector and later on, under the ev_father entry
different chunks have been stored with key image hashes.
Those hashes implements different pHash hashes
( DCT, MW, RD, STRADA).

A vector of string is returned
}
\usage{
VideoSearch(url)
}
\arguments{
\item{url}{URL where the NonSQL database answer queries.}
}
\author{
Silvia Esinosa, <silviaespinosag@gmail.com>,
Joaquin Ordieres, <j.ordieres@upm.es>,
Antonio Bello, <abello@uniovi.es>.
}
\keyword{Video List}
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
# videos<-VideoSearch()
}
