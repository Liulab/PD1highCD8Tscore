ipsmap <-
function (x) {
if (x<=0) {
ips<-0
} else {
if (x>=3) {
 ips<-10
} else {
ips<-round(x*10/3, digits=0)
}
}
return(ips)
}
