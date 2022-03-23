import string
import sys

openfile = open(sys.argv[1], 'r')
newfile = open(sys.argv[2], 'w')

quality_level = 50  # EDIT THIS FOR QUALITY LEVEL DESIRED

num_list = "0123456789"

whole_string = openfile.read()
string_list = whole_string.split("\n");

# newfile.write("Line #\tchr\tnt\tbase\t#reads\tA\tG\tC\tT\tN\ta\tg\tc\tt\tn\t+/-\tI\ti")

line_num = 0

import re

for string in string_list:
	if re.sub("\D", "", string[0:11]) == "":
		digstr = 1
	else:
		digstr = int(re.sub("\D", "", string[0:11]))
	if digstr > 1:
		line_num = line_num + 1
		pre_string = ""
		tab_count = 0
		is_hat = 0
		is_plusminus = 0
		plusminus_val = 0
		plusminus_count = 0
		is_num = 0
		A_ct = 0
		G_ct = 0
		C_ct = 0
		T_ct = 0
		N_ct = 0
		a_ct = 0
		g_ct = 0
		c_ct = 0
		t_ct = 0
		n_ct = 0
		I_ct = 0
		i_ct = 0
		item_array = []
		quality_array = []
		for x in string:
			if x == "\t":
				tab_count = tab_count + 1
			if tab_count < 4:
				pre_string = pre_string + x
			elif tab_count == 4:
				if is_hat:
					is_hat = 0
				elif is_plusminus:
					plusminus_count = plusminus_count + 1
					is_plusminus = 0
					plusminus_val = int(x)
				elif x == "^":
					is_hat = 1
				elif x == "+":
					is_plusminus = 1
				elif x == "-":
					is_plusminus = 1
				elif plusminus_val != 0:
					if plusminus_val == 1:
						if x == "A":
							I_ct = I_ct + 1
						elif x == "G":
							I_ct = I_ct + 1
						elif x == "C":
							I_ct = I_ct + 1
						elif x == "T":
							I_ct = I_ct + 1
						elif x == "N":
							I_ct = I_ct + 1
						elif x == "a":
							i_ct = i_ct + 1
						elif x == "g":
							i_ct = i_ct + 1
						elif x == "c":
							i_ct = i_ct + 1
						elif x == "t":
							i_ct = i_ct + 1
						elif x == "n":
							i_ct = i_ct + 1
						else:
							lalala = 0
						plusminus_val = plusminus_val - 1
					else:
						plusminus_val = plusminus_val - 1
				elif x == "$":
					lalala = 0
				else:
					item_array.append(x)
			elif tab_count == 5:
				quality_array.append(x)
			else:
				lalala = 0
		length_count = 0
		if len(item_array) == len(quality_array):
			for x in item_array:
				item = quality_array[length_count]
				ord_item = ord(item)
				if ord_item >= quality_level:
					if x == "A":
						A_ct = A_ct + 1
					elif x == "G":
						G_ct = G_ct + 1
					elif x == "C":
						C_ct = C_ct + 1
					elif x == "T":
						T_ct = T_ct + 1
					elif x == "N":
						N_ct = N_ct + 1
					elif x == "a":
						a_ct = a_ct + 1
					elif x == "g":
						g_ct = g_ct + 1
					elif x == "c":
						c_ct = c_ct + 1
					elif x == "t":
						t_ct = t_ct + 1
					elif x == "n":
						n_ct = n_ct + 1
					else:
						lalala = 0
				length_count = length_count + 1
		else:
			for x in item_array:
				if x == "A":
					A_ct = A_ct + 1
				elif x == "G":
					G_ct = G_ct + 1
				elif x == "C":
					C_ct = C_ct + 1
				elif x == "T":
					T_ct = T_ct + 1
				elif x == "N":
					N_ct = N_ct + 1
				elif x == "a":
					a_ct = a_ct + 1
				elif x == "g":
					g_ct = g_ct + 1
				elif x == "c":
					c_ct = c_ct + 1
				elif x == "t":
					t_ct = t_ct + 1
				elif x == "n":
					n_ct = n_ct + 1
				else:
					lalala = 0
				length_count = length_count + 1
		newfile.write("\n" + str(line_num) + "\t" + pre_string + "\t" + str(A_ct) + "\t" + str(G_ct) + "\t" + str(
			C_ct) + "\t" + str(T_ct) + "\t" + str(N_ct) + "\t" + str(a_ct) + "\t" + str(g_ct) + "\t" + str(
			c_ct) + "\t" + str(t_ct) + "\t" + str(n_ct) + "\t" + str(plusminus_count) + "\t" + str(I_ct) + "\t" + str(
			i_ct))
	else:
		lalala = 0
newfile.write("\n")




