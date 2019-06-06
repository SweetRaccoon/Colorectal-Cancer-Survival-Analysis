#!/usr/bin/perl -w
use strict;
use warnings;

# This script will open a series of clinical data .xml files downloaded from TCGA and extract
#	clinical data. For an example see
#	/Users/Quincy/Cancer/clinical/gdc_download/gdc_download_20161122_030439/0a7a5131-c59a-4d3e-9ce6-4819499c9116/nationwidechildrens.org_clinical.TCGA-AA-3685.xml

# This script uses a "manifest" file that was downloaded with the clinical data that is 
#	essentially a map of the folder name and file name.

#############################################################################################
# The manifest file should be copied to this directory and named manifest.txt. The file should
#	look like the following

#	id	filename	 md5	 size	state
#	019624ba-47df-4bff-81a7-97feb3088221	nationwidechildrens.org_clinical.TCGA-AG-4021.xml	2ff23be17a86571b1801637bb332c920	31421	live
#	ade620b4-898b-492e-b883-0832652536f7	nationwidechildrens.org_clinical.TCGA-AG-A016.xml	4ba802b6e0cb9d875694512759f4e952	46687	live

# Fill in the location of the folders
my $path = "/Users/Quincy/Cancer/clinical/gdc_download/gdc_download_20161122_030439/";


# The actual XML files for each patient will look like the following
#	<?xml version="1.0" encoding="UTF-8"?>
#	<coad:tcga_bcr xsi:schemaLocation="http://tcga.nci/bcr/xml/clinical/coad/2.7 http://tcga-data.nci.nih.gov/docs/xsd/BCR/tcga.nci/bcr/xml/clinical/coad/2.7/TCGA_BCR.COAD_Clinical.xsd" schemaVersion="2.7" xmlns:coad="http://tcga.nci/bcr/xml/clinical/coad/2.7" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:admin="http://tcga.nci/bcr/xml/administration/2.7" xmlns:clin_shared="http://tcga.nci/bcr/xml/clinical/shared/2.7" xmlns:shared="http://tcga.nci/bcr/xml/shared/2.7" xmlns:shared_stage="http://tcga.nci/bcr/xml/clinical/shared/stage/2.7" xmlns:coad_read_shared="http://tcga.nci/bcr/xml/clinical/shared/coad_read/2.7" xmlns:coad_nte="http://tcga.nci/bcr/xml/clinical/coad/shared/new_tumor_event/2.7/1.0" xmlns:nte="http://tcga.nci/bcr/xml/clinical/shared/new_tumor_event/2.7" xmlns:rx="http://tcga.nci/bcr/xml/clinical/pharmaceutical/2.7" xmlns:rad="http://tcga.nci/bcr/xml/clinical/radiation/2.7" xmlns:follow_up_v1.0="http://tcga.nci/bcr/xml/clinical/coad/followup/2.7/1.0">
#	    <admin:admin>
#	        <admin:bcr xsd_ver="1.17">Nationwide Children's Hospital</admin:bcr>


#############################################################################################
# Open the manifest file and create an arra
open (DATAFILE, "<$path/manifest.txt") || die ("Cannot open file");
my @data = <DATAFILE>;
close DATAFILE;

# Get rid of all newline characters in the sequence array
my $index = 0;
foreach (@data)
	{
	chomp ($data[$index]);
	$index++;
	}

 # get rid of the header line
 shift (@data);
#  foreach (@data) {print "$_\n";}
#############################################################################################

#############################################################################################
# Create matching arrays of folder names and file names and extract patient IDs from file names
my @folder_names; my @file_names; my @patient_ids;

foreach (@data)
	{
	my @temp = split(/\t/,$_);
	push (@folder_names, $temp[0]);
	push (@file_names, $temp[1]);
	}
# foreach (@folder_names) {print "$_\n";}
# foreach (@file_names) {print "$_\n";}


foreach (@file_names)
	{
	my $temp_length = length($_);
	my $temp = substr($_, $temp_length-8, 4);
	push (@patient_ids, $temp);
	}
# foreach (@patient_ids) {print "$_\n";}
#############################################################################################

#############################################################################################
# Open a file a extract a piece of data
my %field_1_hash;
my $field_1_name = "clin_shared:days_to_last_followup";

my %field_2_hash;
my $field_2_name = "clin_shared:days_to_death";

my %field_3_hash;
my $field_3_name = "nte:days_to_new_tumor_event_after_initial_treatment";

my %field_4_hash;
my $field_4_name = "clin_shared:vital_status";

my %field_5_hash;
my $field_5_name = "clin_shared:days_to_last_known_alive";



$index = 0;
foreach (@folder_names)
	{
	my $found = 0;
	# Open each file and create a temporary data array
	my $temp_directory = "$path/$file_names[$index]";
	# print "$temp_directory\n";
	open (TEMPDATA, "<$temp_directory") || die ("Cannot open file");
	my @temp_data = <TEMPDATA>;
	close TEMPDATA;
	
	# Get rid of all newline characters in the sequence array
	my $temp_index = 0;
	foreach (@temp_data)
		{
		chomp ($temp_data[$temp_index]);
		chop ($temp_data[$temp_index]);
		$temp_index++;
		}
	# print "$temp_data[112]\n";
	# foreach (@temp_data) {print "$_\n";}

# Use this loop to temporarily search for something
#	foreach (@temp_data)
#		{
#		if ($_ eq "            <follow_up_v1.0:follow_up version=\"1.0\" sequence=\"3\">")
#			{
#			print "$patient_ids[$index]\n";
#			}
#		}

	
	################
	# field_1
	# This is the subroutine that will find data based on the fields defined above
	$found = 0;
	foreach my $temp_string (@temp_data)
		{
		if ($temp_string =~ /.*<$field_1_name/)
			{
			$found = 1;
			# Test to see if the line has not occurred previously in the file and already been assigned a value
			if (!exists $field_1_hash{$patient_ids[$index]})
				{
				# Test to see if the line ends without a value
				if ($temp_string =~ /\/>$/)
					{
					# If the line ends without a value assign "unknown"
					$field_1_hash{$patient_ids[$index]} = "unknown";
					}
					else
						{					
						# search for the value between the carats
						$temp_string =~ /(>.*<)/;
						# Remove the carats on either side
						my $temp_length = length($&);
						my $temp_trim = substr($&, 1, $temp_length-2);
						# Place the value in the hash
						$field_1_hash{$patient_ids[$index]} = $temp_trim;
						# print "$field_1_hash{$patient_ids[$index]}\n";
						}
				}
				# if it already exists, do the same thing, but add it to the hash using a pipe separator
				else
					{
					# Test to see if the line ends without a value
					if ($temp_string =~ /\/>$/)
						{
						# If the line ends without a value assign "unknown"
						$field_1_hash{$patient_ids[$index]} = "$field_1_hash{$patient_ids[$index]}\|unknown";
						}
						else
							{
							# search for the value between the carats
							$temp_string =~ /(>.*<)/;
							# Remove the carats on either side
							my $temp_length = length($&);
							my $temp_trim = substr($&, 1, $temp_length-2);
							# Place the value in the hash
							$field_1_hash{$patient_ids[$index]} = "$field_1_hash{$patient_ids[$index]}\|$temp_trim";
							}
					}			
			}
		}
	# If the line heading was not found, enter unknown
	if ($found == 0)
		{
		$field_1_hash{$patient_ids[$index]} = "unknown";
		}

	#################


	################
	# field_2
	# This is the subroutine that will find data based on the fields defined above
	$found = 0;
	foreach my $temp_string (@temp_data)
		{
		if ($temp_string =~ /.*<$field_2_name/)
			{
			$found = 1;
			# Test to see if the line has not occurred previously in the file and already been assigned a value
			if (!exists $field_2_hash{$patient_ids[$index]})
				{
				# Test to see if the line ends without a value
				if ($temp_string =~ /\/>$/)
					{
					# If the line ends without a value assign "unknown"
					$field_2_hash{$patient_ids[$index]} = "unknown";
					}
					else
						{					
						# search for the value between the carats
						$temp_string =~ /(>.*<)/;
						# Remove the carats on either side
						my $temp_length = length($&);
						my $temp_trim = substr($&, 1, $temp_length-2);
						# Place the value in the hash
						$field_2_hash{$patient_ids[$index]} = $temp_trim;
						# print "$field_2_hash{$patient_ids[$index]}\n";
						}
				}
				# if it already exists, do the same thing, but add it to the hash using a pipe separator
				else
					{
					# Test to see if the line ends without a value
					if ($temp_string =~ /\/>$/)
						{
						# If the line ends without a value assign "unknown"
						$field_2_hash{$patient_ids[$index]} = "$field_2_hash{$patient_ids[$index]}\|unknown";
						}
						else
							{
							# search for the value between the carats
							$temp_string =~ /(>.*<)/;
							# Remove the carats on either side
							my $temp_length = length($&);
							my $temp_trim = substr($&, 1, $temp_length-2);
							# Place the value in the hash
							$field_2_hash{$patient_ids[$index]} = "$field_2_hash{$patient_ids[$index]}\|$temp_trim";
							}
					}			
			}
		}
	# If the line heading was not found, enter unknown
	if ($found == 0)
		{
		$field_2_hash{$patient_ids[$index]} = "unknown";
		}

	#################
	
	################
	# field_3
	# This is the subroutine that will find data based on the fields defined above
	$found = 0;
	foreach my $temp_string (@temp_data)
		{
		if ($temp_string =~ /.*<$field_3_name/)
			{
			$found = 1;
			# Test to see if the line has not occurred previously in the file and already been assigned a value
			if (!exists $field_3_hash{$patient_ids[$index]})
				{
				# Test to see if the line ends without a value
				if ($temp_string =~ /\/>$/)
					{
					# If the line ends without a value assign "unknown"
					$field_3_hash{$patient_ids[$index]} = "unknown";
					}
					else
						{					
						# search for the value between the carats
						$temp_string =~ /(>.*<)/;
						# Remove the carats on either side
						my $temp_length = length($&);
						my $temp_trim = substr($&, 1, $temp_length-2);
						# Place the value in the hash
						$field_3_hash{$patient_ids[$index]} = $temp_trim;
						# print "$field_3_hash{$patient_ids[$index]}\n";
						}
				}
				# if it already exists, do the same thing, but add it to the hash using a pipe separator
				else
					{
					# Test to see if the line ends without a value
					if ($temp_string =~ /\/>$/)
						{
						# If the line ends without a value assign "unknown"
						$field_3_hash{$patient_ids[$index]} = "$field_3_hash{$patient_ids[$index]}\|unknown";
						}
						else
							{
							# search for the value between the carats
							$temp_string =~ /(>.*<)/;
							# Remove the carats on either side
							my $temp_length = length($&);
							my $temp_trim = substr($&, 1, $temp_length-2);
							# Place the value in the hash
							$field_3_hash{$patient_ids[$index]} = "$field_3_hash{$patient_ids[$index]}\|$temp_trim";
							}
					}			
			}
		}
	# If the line heading was not found, enter unknown
	if ($found == 0)
		{
		$field_3_hash{$patient_ids[$index]} = "unknown";
		}

	#################
	
	################
	# field_4
	# This is the subroutine that will find data based on the fields defined above
	$found = 0;
	foreach my $temp_string (@temp_data)
		{
		if ($temp_string =~ /.*<$field_4_name/)
			{
			$found = 1;
			# Test to see if the line has not occurred previously in the file and already been assigned a value
			if (!exists $field_4_hash{$patient_ids[$index]})
				{
				# Test to see if the line ends without a value
				if ($temp_string =~ /\/>$/)
					{
					# If the line ends without a value assign "unknown"
					$field_4_hash{$patient_ids[$index]} = "unknown";
					}
					else
						{					
						# search for the value between the carats
						$temp_string =~ /(>.*<)/;
						# Remove the carats on either side
						my $temp_length = length($&);
						my $temp_trim = substr($&, 1, $temp_length-2);
						# Place the value in the hash
						$field_4_hash{$patient_ids[$index]} = $temp_trim;
						# print "$field_4_hash{$patient_ids[$index]}\n";
						}
				}
				# if it already exists, do the same thing, but add it to the hash using a pipe separator
				else
					{
					# Test to see if the line ends without a value
					if ($temp_string =~ /\/>$/)
						{
						# If the line ends without a value assign "unknown"
						$field_4_hash{$patient_ids[$index]} = "$field_4_hash{$patient_ids[$index]}\|unknown";
						}
						else
							{
							# search for the value between the carats
							$temp_string =~ /(>.*<)/;
							# Remove the carats on either side
							my $temp_length = length($&);
							my $temp_trim = substr($&, 1, $temp_length-2);
							# Place the value in the hash
							$field_4_hash{$patient_ids[$index]} = "$field_4_hash{$patient_ids[$index]}\|$temp_trim";
							}
					}			
			}
		}
	# If the line heading was not found, enter unknown
	if ($found == 0)
		{
		$field_4_hash{$patient_ids[$index]} = "unknown";
		}

	#################
	
	################
	# field_5
	# This is the subroutine that will find data based on the fields defined above
	$found = 0;
	foreach my $temp_string (@temp_data)
		{
		if ($temp_string =~ /.*<$field_5_name/)
			{
			$found = 1;
			# Test to see if the line has not occurred previously in the file and already been assigned a value
			if (!exists $field_5_hash{$patient_ids[$index]})
				{
				# Test to see if the line ends without a value
				if ($temp_string =~ /\/>$/)
					{
					# If the line ends without a value assign "unknown"
					$field_5_hash{$patient_ids[$index]} = "unknown";
					}
					else
						{					
						# search for the value between the carats
						$temp_string =~ /(>.*<)/;
						# Remove the carats on either side
						my $temp_length = length($&);
						my $temp_trim = substr($&, 1, $temp_length-2);
						# Place the value in the hash
						$field_5_hash{$patient_ids[$index]} = $temp_trim;
						# print "$field_5_hash{$patient_ids[$index]}\n";
						}
				}
				# if it already exists, do the same thing, but add it to the hash using a pipe separator
				else
					{
					# Test to see if the line ends without a value
					if ($temp_string =~ /\/>$/)
						{
						# If the line ends without a value assign "unknown"
						$field_5_hash{$patient_ids[$index]} = "$field_5_hash{$patient_ids[$index]}\|unknown";
						}
						else
							{
							# search for the value between the carats
							$temp_string =~ /(>.*<)/;
							# Remove the carats on either side
							my $temp_length = length($&);
							my $temp_trim = substr($&, 1, $temp_length-2);
							# Place the value in the hash
							$field_5_hash{$patient_ids[$index]} = "$field_5_hash{$patient_ids[$index]}\|$temp_trim";
							}
					}			
			}
		}
	# If the line heading was not found, enter unknown
	if ($found == 0)
		{
		$field_5_hash{$patient_ids[$index]} = "unknown";
		}

	#################

	$index++;
	}
# print map { "$_\t$field_1_hash{$_}\n" } keys %field_1_hash;

print "patient_id\t$field_1_name\t$field_2_name\t$field_3_name\t$field_4_name\t$field_5_name\n";
foreach (@patient_ids)
	{
	print "$_\t$field_1_hash{$_}\t$field_2_hash{$_}\t$field_3_hash{$_}\t$field_4_hash{$_}\t$field_5_hash{$_}\n";
	}



