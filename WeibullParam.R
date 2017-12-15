#---------------------------------Estimación Parámetros de una Distribución Weibull----------------------------------------------------------------

#V1 es el vector de velocidades medido

#k son los valores que se iterarán

#method puede ser: LS (Least Squares), Mom (Moments), Wasp, MLE (Maximum Likelihood), MMLE (Modified Maximum Likelihood)

#Sólo puede utilizarse un método cada vez que se llama a la función!

WeibullParam<-function(V1,k=seq(1,2,0.01), method){
        
        library(ggthemes)
        library(extrafont)
        library(ggplot2)
        library(gridExtra)
        
        barfill <- "#4271AE"
        barlines <- "#1F3552"
        
        V1<-V1[V1>0]
        
        V3<-floor(V1)
        
        if(is.vector(V1)==FALSE){
                stop("V1 debe ser un vector")
        }
        
        if(length(which(V1>90))!=0){
                stop("V1 debe ser filtrado. Valores superiores a 90 encontrados")
        }
        
        if(is.numeric(V1)==FALSE){
                stop("V1 debe ser numérico")
        }
        
        if(is.character(method)==FALSE){
                stop("method es un argumento inválido")
        }
        
        if(is.vector(k)==FALSE | is.numeric(k)==FALSE){
                stop("k debe ser un vector numérico")
        }
        
        if(method=="LS"){
                
                V1<-V1[V1>0 & is.na(V1)==FALSE]
                
                ordenar<-order(V1)
                
                V1<-V1[ordenar]
                frecuencia<-seq(1,length(V1),1)/(length(V1)+1)
                
                y<-log(-log(1-frecuencia))
                x<-log(V1)
                
                rx<-lm(y~x)
                
                df1<-as.data.frame(matrix(nrow = length(x), ncol = 2)); df1$x<-x; df1$y<-y
                
                p6<-ggplot(data = df1, aes(x=x,y=y))+
                        geom_point(size=3, colour="red", fill=I("red"))+
                        geom_smooth(method = lm, formula = y~x, lwd=1.5)+
                        scale_x_continuous(name = "ln(U)",
                                           limits=c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))) +
                        scale_y_continuous(name = "ln(ln(1/1-F(u)))", limits = c(min(y, na.rm = TRUE),max(y, na.rm = TRUE))) +
                        labs(colour="Ajuste")+
                        ggtitle("Estimación de parámetros por método gráfico") +
                        theme_economist() +
                        theme(legend.position = "bottom", legend.direction = "horizontal",
                              legend.box = "horizontal",
                              legend.key.size = unit(1, "cm"),
                              plot.title = element_text(family="Bauhaus 93", size = 18),
                              text = element_text(family = "Bauhaus 93"),
                              axis.title = element_text(size = 12, face = "bold"),
                              legend.text = element_text(size = 9),
                              legend.title=element_text(face = "bold", size = 9, family = "Bauhaus 93"))
                
                k<-as.vector(rx$coefficients[2])
                intercept<-as.vector(rx$coefficients[1])
                A<-exp(-intercept/k)
                
                #Cálculo de Chi_Cuadrado y RMSE:
                
                V1<-V1[is.na(V1)==FALSE]; V3<-floor(V1)
                
                Vmedia1<-mean(V1)
                V4<-as.vector(table(V3))
                suma5<-sum(V4)
                frecuencia<-V4/suma5
                Vbinmedio<-as.vector(tapply(V1,V3,mean))
                
                frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
                
                for(i in 1:length(frecuencia)){
                        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
                }
                
                suma_chi<-0
                
                for(i in 1:length(frecuencia1)){
                        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
                }
                
                media<-mean(frecuencia, na.rm = TRUE)
                suma11<-0
                
                for(i in 1:length(frecuencia1)){
                        suma11<-suma11+(frecuencia[i]-media)^2
                }
                
                R<-1-(suma_chi/suma11)
                
                Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
                
                RMSE<-(suma_chi/length(Vbinmedio))^0.5
                
                
                Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE,R_2_Weibull=R, R_2_Ajuste=summary(rx)$r.squared)
                return(Coeficientes)
        }
        
        if(method=="Mom"){
                
                vmedia<-mean(V1[is.na(V1)==FALSE])
                smedia<-sd(V1[is.na(V1)==FALSE])
                
                numerador<-vector(mode = "numeric", length = length(k))
                denominador<-vector(mode="numeric", length = length(k))
                formula<-vector(mode = "numeric", length = length(k))
                
                for(j in 1:length(k)){
                        numerador[j]<-(smedia^2)*(gamma(1+(1/k[j]))^2)
                        denominador[j]<-(vmedia^2)*(gamma(1+(2/k[j]))-gamma(1+(1/k[j]))^2)
                        formula[j]<-(numerador[j]/denominador[j])-1
                }
                
                minimo<-min(abs(formula))
                loc_min<-which(abs(formula)==minimo)
                
                k<-k[loc_min]
                A<-vmedia/gamma(1+(1/k))
                
                #Cálculo de Chi_Cuadrado y RMSE:
                
                Vmedia1<-mean(V1[is.na(V1)==FALSE])
                V4<-as.vector(table(V3))
                suma5<-sum(V4)
                frecuencia<-V4/suma5
                Vbinmedio<-as.vector(tapply(V1,V3,mean))
                
                frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
                
                for(i in 1:length(frecuencia)){
                        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
                }
                
                suma_chi<-0
                
                for(i in 1:length(frecuencia1)){
                        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
                }
                
                media<-mean(frecuencia, na.rm = TRUE)
                suma11<-0
                
                for(i in 1:length(frecuencia1)){
                        suma11<-suma11+(frecuencia[i]-media)^2
                }
                
                R<-1-(suma_chi/suma11)
                
                Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
                
                RMSE<-(suma_chi/length(Vbinmedio))^0.5
                
                
                
                z<-as.numeric(names(tapply(V1,V3,mean)))
                
                frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")
                
                frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; frecuencia1_df$frecuencia1<-frecuencia1
                
                
                
                par(mfrow=c(1,1))
                
                barfill <- "#4271AE"
                barlines <- "#1F3552"
                
                p7 <- ggplot(frecuencia1_df) +
                        geom_col(aes(x =z+0.5, y=frecuencia),
                                 colour = barlines, fill = barfill) +
                        scale_x_continuous(name = "Velocidad (m/s)",
                                           breaks = seq(1,length(frecuencia1),1),
                                           limits=c(0, length(frecuencia))) +
                        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
                        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
                        labs(colour="Ajuste")+
                        ggtitle("Distribución ajustada a Weibull por Momentos") +
                        theme_economist() +
                        theme(legend.position = c(0.75,0.75), legend.direction = "horizontal",
                              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
                              legend.box = "horizontal",
                              legend.key.size = unit(1, "cm"),
                              plot.title = element_text(family="comic-sans", size = 22),
                              text = element_text(family = "comic-sans"),
                              axis.title = element_text(size = 12, face = "bold"),
                              legend.text = element_text(size = 9),
                              legend.title=element_blank())+
                        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))
                
                p7
                
                Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p7)
                
                
                
                return(Coeficientes)
        }
        
        if(method=="Wasp"){
                
                V1<-V1[is.na(V1)==FALSE]
                
                ordenar<-order(V1); V1<-V1[ordenar]; V3<-floor(V1)
                
                V2<-abs(V1-mean(V1, na.rm = TRUE))
                
                loc<-which(V2==min(V2, na.rm = TRUE))
                
                X<-1-(loc[length(loc)]/(length(V1)+1))
                
                Vmedia<-mean(V1, na.rm = TRUE); suma<-sum(V1^3, na.rm = TRUE); N<-length(V1)
                
                N<-length(V1[is.na(V1)==FALSE])
                Vmedia<-mean(V1[is.na(V1)==FALSE])
                
                suma<-0
                
                for(i in 1:N){
                        if(is.na(V1[i])==FALSE){
                                
                                suma<-suma+(V1[i])^3
                        }
                }
                
                #Cálculo de frecuencias:
                V4<-as.vector(table(V3))
                suma1<-sum(V4)
                frecuencia<-V4/suma1
                
                 #Aplicamos la fórmula
                
                numerador1<-vector(mode = "numeric", length = length(k))
                numerador2<-vector(mode = "numeric", length = length(k))
                formula<-vector(mode = "numeric", length = length(k))
                
                for(j in 1:length(k)){
                        numerador1[j]<-(-log(X))^(1/k[j])
                        numerador2[j]<-(suma/(N*gamma(1+(3/k[j]))))^(1/3)
                        
                        formula[j]<-((numerador1[j]*numerador2[j])/Vmedia)-1
                }
                
                #Buscamos la k para la que el resultado está más próximo a 0:
                
                minimo<-min(abs(formula))
                loc_min<-which(abs(formula)==minimo)
                k<-k[loc_min]
                
                
                A<-(suma/(N*gamma(1+(3/k))))^(1/3)
                
                #Cálculo de Chi_Cuadrado y RMSE:
                
                Vbinmedio<-as.vector(tapply(V1, V3, mean))
                frecuencia1<-vector(mode = "numeric", length = length(Vbinmedio))
                
                for(i in 1:length(Vbinmedio)){
                        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
                }
                #
                suma_chi<-0
                #
                for(i in 1:length(frecuencia1)){
                        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
                }
                
                media<-mean(frecuencia, na.rm = TRUE)
                suma11<-0
                #
                for(i in 1:length(frecuencia1)){
                        suma11<-suma11+(frecuencia[i]-media)^2
                }
                #
                R<-1-(suma_chi/suma11)
                #
                Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
                #
                RMSE<-(suma_chi/length(Vbinmedio))^0.5
                
                z<-as.numeric(names(tapply(V1,V3,mean)))
                
                par(mfrow=c(1,1))
                
                frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")
                
                frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; frecuencia1_df$frecuencia1<-frecuencia1
                
                barfill <- "#4271AE"
                barlines <- "#1F3552"
                
                p7 <- ggplot(frecuencia1_df) +
                        geom_col(aes(x =z+0.5, y=frecuencia),
                                 colour = barlines, fill = barfill) +
                        scale_x_continuous(name = "Velocidad (m/s)",
                                           breaks = seq(1,length(frecuencia1),1),
                                           limits=c(0, length(frecuencia))) +
                        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
                        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
                        labs(colour="Ajuste")+
                        ggtitle("Distribución ajustada a Weibull por Wasp") +
                        theme_economist() +
                        theme(legend.position = c(0.75, 0.75), legend.direction = "horizontal",
                              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
                              legend.box = "horizontal",
                              legend.key.size = unit(1, "cm"),
                              plot.title = element_text(family="comic-sans", size = 22),
                              text = element_text(family = "comic-sans"),
                              axis.title = element_text(size = 12, face = "bold"),
                              legend.text = element_text(size = 9),
                              legend.title=element_blank())+
                        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))
                
                p7
                
                Coeficientes<-list(A=A, k=k, RMSE=RMSE, Chi_Cuadrado=Chi_Cuadrado, R2=R, graph=p7)
                
                
                return(Coeficientes)
        }
        
        if(method=="MLE"){
                
                N<-length(V1[is.na(V1)==FALSE])
                
                suma<-0
                
                for(i in 1:N){
                        if(is.na(V1[i])==FALSE){
                                suma<-suma+log(V1[i])
                        }
                }
                
                suma1<-0
                
                for(i in 1:N){
                        if(is.na(V1[i])==FALSE){
                                suma1<-suma1+(V1[i])^k
                        }
                }
                
                suma2<-0
                
                for(i in 1:N){
                        if(is.na(V1[i])==FALSE){
                                suma2<-suma2+(((V1[i])^k)*log(V1[i]))
                        }
                }
                
                #Calculamos la función:
                numerador1<-k*N*suma2
                numerador2<-k*suma1*suma
                denominador<-N*suma1
                formula<-((numerador1-numerador2)/denominador)-1
                
                #Obtenemos k en base al valor más próximo a 0 obtenido en la fórmula:
                minimo<-min(abs(formula))
                loc_min<-which(abs(formula)==minimo)
                k<-k[loc_min]
                
                
                #Obtenemos A a partir de k:
                suma4<-0
                
                for(i in 1:N){
                        if(is.na(V1[i])==FALSE){
                                suma4<-suma4+((V1[i])^k)
                        }
                }
                
                A<-(suma4/N)^(1/k)
                
                #Cálculo de Chi_Cuadrado y RMSE:
                
                Vmedia<-mean(V1[is.na(V1)==FALSE])
                V4<-as.vector(table(V3))
                suma5<-sum(V4)
                frecuencia<-V4/suma5
                Vbinmedio<-as.vector(tapply(V1,V3,mean))
                
                frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
                
                for(i in 1:length(frecuencia)){
                        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
                }
                
                suma_chi<-0
                
                for(i in 1:length(frecuencia1)){
                        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
                }
                
                media<-mean(frecuencia, na.rm = TRUE)
                suma11<-0
                
                for(i in 1:length(frecuencia1)){
                        suma11<-suma11+(frecuencia[i]-media)^2
                }
                
                R<-1-(suma_chi/suma11)
                
                Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
                
                RMSE<-(suma_chi/length(Vbinmedio))^0.5
                
                z<-as.numeric(names(tapply(V1,V3,mean)))
                
                par(mfrow=c(1,1))
                
                frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")
                
                frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; frecuencia1_df$frecuencia1<-frecuencia1
                
                p7 <- ggplot(frecuencia1_df) +
                        geom_col(aes(x =z+0.5, y=frecuencia),
                                 colour = barlines, fill = barfill) +
                        scale_x_continuous(name = "Velocidad (m/s)",
                                           breaks = seq(1,length(frecuencia1),1),
                                           limits=c(0, length(frecuencia))) +
                        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
                        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
                        labs(colour="Ajuste")+
                        ggtitle("Distribución ajustada a Weibull por Máxima Verosimilitud") +
                        theme_economist() +
                        theme(legend.position = c(0.75, 0.75), legend.direction = "horizontal",
                              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
                              legend.box = "horizontal",
                              legend.key.size = unit(1, "cm"),
                              plot.title = element_text(family="Bauhaus 93", size = 18),
                              text = element_text(family = "comic-sans"),
                              axis.title = element_text(size = 12, face = "bold"),
                              legend.text = element_text(size = 9),
                              legend.title=element_blank())+
                        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))
                
                p7
                
                Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p7)
                
                return(Coeficientes)
        }
        
        if(method=="MMLE"){
                
                N<-length(V1[is.na(V1)==FALSE])
                
                suma<-0
                
                for(i in 1:N) {
                        if(is.na(V1[i])==FALSE){
                                suma<-suma+(log(V1[i])*log(V1[i]))
                                
                        }
                }
                
                suma1<-0
                
                for(i in 1:N){
                        if(is.na(V1[i])==FALSE){
                                suma1<-suma1+log(V1[i])
                        }
                }
                
                numerador<-N*(N-1)
                denominador<-(N*suma)-((suma1)^2)
                k<-(pi/sqrt(6))*(numerador/denominador)^0.5
                
                suma3<-0
                
                for(i in 1:N){
                        if(is.na(V1[i])==FALSE){
                                suma3<-suma3+(V1[i])^k
                        }
                }
                
                A<-(suma3/N)^(1/k)
                
                #Cálculo de Chi_Cuadrado y RMSE:
                
                Vmedia1<-mean(V1[is.na(V1)==FALSE])
                V4<-as.vector(table(V3))
                suma5<-sum(V4)
                frecuencia<-V4/suma5
                Vbinmedio<-as.vector(tapply(V1,V3,mean))
                
                frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
                
                for(i in 1:length(frecuencia)){
                        frecuencia1[i]<-(k/A)*((Vbinmedio[i]/A)^(k-1))*exp(-(Vbinmedio[i]/A)^k)
                }
                
                suma_chi<-0
                
                for(i in 1:length(frecuencia1)){
                        suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
                }
                
                media<-mean(frecuencia, na.rm = TRUE)
                suma11<-0
                
                for(i in 1:length(frecuencia1)){
                        suma11<-suma11+(frecuencia[i]-media)^2
                }
                
                R<-1-(suma_chi/suma11)
                
                Chi_Cuadrado<-suma_chi/(length(Vbinmedio)-2)
                
                RMSE<-(suma_chi/length(Vbinmedio))^0.5
                
                z<-as.numeric(names(tapply(V1,V3,mean)))
                
                par(mfrow=c(1,1))
                
                frecuencia1_df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3)); colnames(frecuencia1_df)<-c("frecuencia", "z", "frecuencia1")
                
                frecuencia1_df$frecuencia<-frecuencia; frecuencia1_df$z<-z; frecuencia1_df$frecuencia1<-frecuencia1
                
                p7 <- ggplot(frecuencia1_df) +
                        geom_col(aes(x =z+0.5, y=frecuencia),
                                 colour = barlines, fill = barfill) +
                        scale_x_continuous(name = "Velocidad (m/s)",
                                           breaks = seq(1,length(frecuencia1),1),
                                           limits=c(0, length(frecuencia))) +
                        scale_y_continuous(name = "Frecuencia", limits = c(0,0.15)) +
                        geom_line(data = frecuencia1_df, aes(x=z+0.5, y=frecuencia1, colour="red"), size=1.5)+
                        labs(colour="Ajuste")+
                        ggtitle("Distribución ajustada a Weibull por Modificación Máxima Verosimilitud") +
                        theme_economist() +
                        theme(legend.position = c(0.75, 0.75), legend.direction = "horizontal",
                              legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
                              legend.box = "horizontal",
                              legend.key.size = unit(1, "cm"),
                              plot.title = element_text(family="comic-sans", size = 16),
                              text = element_text(family = "comic-sans"),
                              axis.title = element_text(size = 12, face = "bold"),
                              legend.text = element_text(size = 9),
                              legend.title=element_blank())+
                        scale_colour_discrete("Point",labels=c("Ajuste Weibull"))
                
                p7
                

                Coeficientes<-list(A=A, k=k, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p7)
                
                
                return(Coeficientes)
        }
        
        else{
                stop("method es un argumento inválido")
        }
}
