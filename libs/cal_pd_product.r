cal_pd_product <- function(xs, h, ft, sc = 1, error = 1, mask = NULL, FUN = function(x) x) {
    ht = FUN(h)
    yay = ft
    if (is.null(mask)) browser()
    ftt = FUN((ft+h*sc)/sc)
    #ftt = addLayer(ftt - error*1.5, ftt + error*1.5)
    FUN4cell <- function(x, y) {
        p = dnorm(xs, x, error*0.2)
        q = dnorm(xs, y, error*0.2)
        #browser()
        mean(sqrt(p*q))/mean(p)#log(p/q))
    }
    out = h
    pv = pv0 = mapply(FUN4cell, ht[mask], ftt[mask])
    
    pv = 1-pv#exp(-pv)
    out[mask] = pv
    test = ft < 0
    out[test] = -out[test]   
    #browser()
    out
}    
