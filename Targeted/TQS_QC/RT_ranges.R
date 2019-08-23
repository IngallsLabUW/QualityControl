df = data.frame(x=c(2,-1,5,4,7,8), y=c(30521, 1800, 25000,1000000, -5, 10))
limits = data.frame("var"=c("x", "y"), min=c(1,0), max=c(5,99999))

sweep(df, 2, limits[, 2], FUN='>') & sweep(df, 2, limits[, 3], FUN='<')


df = RT.flags.added
limits = data.frame(sample = RT.flags.added$Replicate.Name, min=RT.flags.added$RT.min, max=RT.flags.added$RT.max)
     

#
sweep(df, 10, limits[, 2], FUN='>', if(NA) {}) & sweep(df, 10, limits[, 3], FUN='<', if (NA) {})


#
subset(df, df$Retention.Time > (df$RT.min - RT.flex) & df$Retention.Time < (df$RT.max + RT.flex))
