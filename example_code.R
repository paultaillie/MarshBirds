mini_iris <- iris[c(1, 51, 101), ]

gather(mini_iris, key = flower_att, value = measurement,
       Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)



df <- 
  structure(list(id = c(2L, 2L, 4L, 5L, 5L, 5L, 5L), start_end = structure(c(2L, 
                                                                             1L, 2L, 2L, 1L, 2L, 1L), .Label = c("end", "start"), class = "factor"), 
                 date = structure(c(6L, 7L, 3L, 8L, 9L, 10L, 11L), .Label = c("1979-01-03", 
                                                                              "1979-06-21", "1979-07-18", "1989-09-12", "1991-01-04", "1994-05-01", 
                                                                              "1996-11-04", "2005-02-01", "2009-09-17", "2010-10-01", "2012-10-06"
                 ), class = "factor")), .Names = c("id", "start_end", "date"
                 ), row.names = c(3L, 4L, 7L, 8L, 9L, 10L, 11L), class = "data.frame")

df %>%
  group_by(start_end, id) %>%
  mutate(ind = row_number()) %>%
  spread(start_end, date) %>% 
  select(start, end)



