import numpy as np
import streamlit as st
#from Codon_opt_functions import *
from io import StringIO
import re
import pandas as pd
from codon_opt_functions import *
import numpy as np
from PIL import Image
import base64
import time



st.title('Ultimate Codon Optimizer >9000')
image = Image.open('vegeta.jpg')
st.image(image, caption = 'For real')
#file_ = open("its-over9000-vegeta.gif", "rb")
#contents = file_.read()
#data_url = base64.b64encode(contents).decode("utf-8")
#file_.close()

#st.markdown(
#    f'<img src="data:image/gif;base64,{data_url}" alt="cat gif">',
#    unsafe_allow_html=True,
#)
st.text('This is a web app which allows the codon optimization of a \ngiven Amino acid sequence and codon usage table')



st.header('Input:')
# Upload Profile
profile_file = st.file_uploader('Codon Usage Profile')

if profile_file is not None:

	profile_stringio = StringIO(profile_file.getvalue().decode("utf-8"))
	# To read file as string:
	profile_string = profile_stringio.read()
	profile = parse(profile_string)

#st.header('Amino acid sequence Input:')
AA_seq = st.text_input("Amino acid Sequnce:")

if len(AA_seq) != 0:
	if  profile_file is not None:
		#st.write(type(AA_seq), len(AA_seq))
		codon_opt_df = codon_opt(AA_seq = AA_seq, profile = profile)
	else:
		codon_opt_df = None
else:
	codon_opt_df = None
# Create a section for Profile header
st.header('Optimized Sequence')
if codon_opt_df is not None:
	total_mean_dist = mean_dist_fitness_func_arr(np.asarray(codon_opt_df))
	st.write("The toal mean distance is %s" %(total_mean_dist))
	st.write(codon_opt_df)
	csv = codon_opt_df.to_csv(index=False).encode('utf-8')
	st.download_button(
		"Press to Download",
		csv,
		"file.csv",
		"text/csv",
		key='download-csv'
	)

#input = st.text_area("Give Your Input","",height =20)
#modelresponse = model_function(input)
#st.text_area(label ="",value=input, height =100)
	#st.download_button('Download CSV', codon_opt_df, 'text/csv')

st.markdown("""
<style>
.stTextArea [data-baseweb=base-input] {
font-size:120%; 
font-weight:bold; 
-webkit-text-fill-color: #22fb4f; 
background-color: black;
border: 2px;
border-radius: 3px;
} 
</style>
""", unsafe_allow_html=True)

output_empty = ""
output_1 = "...\n"
output_2 = "> establishhing VPN Tunnel to Москва, Россия\n"
output_3 = "> establishhing VPN Tunnel to Москва, Россия...\n"
output_4 = "> establishhing VPN Tunnel to Москва, Россия"
output_5 = "> Connection secured\n"
output_6 = "> Client: FSB Vladi is logged in now\n"
output_7 = "> Payload delivered\n"
output_8 = "> Starting encryption"



time.sleep(0.5)

st.text_area(
    label="Text area:",
    value=output_1 + output_2 + output_3,
    height=300)