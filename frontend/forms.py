from django import forms
DATABASE_CHOICES = (
	("nucleotide","Nucleotide"),
	("protein","Protein"),
	("pubmed","Pubmed"),
	)

class SubmitForm(forms.Form):
	email = forms.EmailField()
	search_term = forms.CharField()
	species_file = forms.FileField()
	database = forms.ChoiceField(choices=DATABASE_CHOICES)

