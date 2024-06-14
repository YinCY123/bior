####HSV####
# 计算色调  
calculate_hue <- function(r, g, b) {  
  # 归一化RGB值  
  r <- r / 255  
  g <- g / 255  
  b <- b / 255  
  
  # 确定RGB中的最大值、中间值和最小值  
  max.val <- max(r, g, b)  
  min.val <- min(r, g, b)  
  diff <- max.val - min.val  
  
  # 计算色调  
  if (max.val == min.val) {  
    hue <- 0 # 灰色，无色调  
  } else if (max.val == r) {  
    hue <- (60 * ((g - b) / diff) + 360) %% 360  
  } else if (max.val == g) {  
    hue <- (60 * ((b - r) / diff) + 120) %% 360  
  } else if (max.val == b) {  
    hue <- (60 * ((r - g) / diff) + 240) %% 360  
  }  
  
  return(hue)  
}  

#确定色调所属范围
determine_color <- function(hue) {  
  if (hue < 0 || hue > 360) {  
    stop("Hue value must be between 0 and 360.")  
  }  
  
  if(hue==0){
    return("gray")
  } else if (hue > 0 && hue < 10 || hue >= 345 && hue < 360) {  
    return("Red")  
  } else if (hue >= 10 && hue < 25) {  
    return("brown")  
  } else if (hue >= 25 && hue < 50) {  
    return("Orange")  
  } else if (hue >= 50 && hue < 75) {  
    return("Yellow")  
  } else if (hue >= 75 && hue < 105) {  
    return("Yellow Green") #  
  } else if (hue >= 105 && hue < 135) {  
    return("Green")  
  } else if (hue >= 135 && hue < 140) {  
    return("Cyan") # 青绿色  
  } else if (hue >= 140 && hue < 175) {  
    return("Green") 
  } else if (hue >= 175 && hue < 255) {  
    return("Blue")  
  } else if (hue >= 255 && hue < 285) {  
    return("Purple")  
  } else if (hue >= 285 && hue < 315) {  
    return("purple") 
  } else if (hue >= 315 && hue < 345) {  
    return("Pink")  
  }  
}


color_to_name<-function(rgb_value){
  color_chinese=data.frame(row.names = c("gray","Red","brown","Orange","Yellow","Yellow Green","Green","Cyan","Blue","Purple","Pink"),
                           name=c("灰色","红色","棕色","橙色","黄色","黄绿色","绿色","青色","蓝色","紫色","粉红色"))
  if (is.character(rgb_value)) { rgb_value <- col2rgb(rgb_value) }
  hue <- calculate_hue(rgb_value[1],rgb_value[2],rgb_value[3])
  color_en<-determine_color(hue)
  color_cn<-color_chinese[color_en,1]
  return(c(color_en,color_cn))
}