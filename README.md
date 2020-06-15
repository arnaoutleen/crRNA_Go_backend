# Overview
crRNA Go is a a machine-learning-based web application that can help predict the location of energetically-favorable binding sites for Cas13b on an ssRNA viral genome, and associated energies of such locations. Such ML application should be able to be trained by the user. The web-app will use these sequences to design a crRNA template that can be ordered and used in the lab. 

This repository has the back-end code that includes the code to generate test data, the code for the different types of ML algorithms that I tried out, and the code to test out the ML algorithm, along with a sample bidirectional LSTM model I made and a test dataset (.csv)

For the full web-app, please go to my other repo: https://github.com/arnaoutleen/crRNA_Go


# To better understand how this program works
https://docs.google.com/presentation/d/10fOaQTVwgwUEoDLApQ1MJ4ZB8ncaO7QaIXrIc_JAOlc/edit?usp=sharing


# Contents

|

|-```README.md```: obvious what it does

|-```load_test_data.txt```: link to where test_data.csv is online, to download it

|-```lstm_bidirectional.py```: bidirectional LSTM

|-```lstm.py```: regular LSTM

|-```genome_draw.py```: making graph of where binding sites are on the genome of virus

|-```backend_load_predict_csvwrite_draw.py```: test for flask app, combining all functiona together

|-```print_crRNA_string.py```: prototype of chassis thing

|-```test_data_generator.py```: generates your local copy of test_data.csv, new test data

|-```trying_new_lstm_params.py```: made a second copy of lstm_bidirectional.py

|-```try_one_hotting.py```: one-hot encoder

|-```my_model_bidirectional.h5```: sample ML model

|-```Test Files```:

    |-```load_test_data.txt```: link to where test_data.csv is online, to download it
  
    |-``` my_model_bidirectional.h5```: sample ML model
  

# Co-Dependencies
 to run, this project requires: Python 3.7.6, Flask 1.1.1, Werkzeug 1.0.0, TensorFlow 2.2.0, Keras 2.3.1, Pandas 1.0.1, PIL/pillow 7.1.2, numpy 1.18.1, Biopython 1.77 (this includes Bio.Seq and Bio.SeqFeature and Bio.Graphics), sklearn 0.0, time

 It was built and run on Ubuntu 18.04


# How to run
Download all the files into a root directory. The root directory should include all the .py files along with a folder for "Test Files". Download the test_data.csv files as instructed in "load_test_data.txt" and add it to Test Files and the root directory. (Alternatively you can run ``` test_data_generator.py``` and make your own local different test_data.csv) Navigate your terminal window to the root directory, then to run any python script, type in terminal:

```
python (scriptname).py
```


# License
I mean, this was an app developed for a college class so ... yeah use it in whatever way you'd like just please credit me and link to this repo so that I can see what people do with this.
