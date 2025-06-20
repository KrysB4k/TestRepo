
function resetText(ttip)
{
  ttip.innerHTML = "Copy to clipboard";
}

function copyCode(btn)
{
	let tooltip = btn.firstElementChild;
	tooltip.innerHTML = "Text copied!";
	setTimeout(resetText, 1500, tooltip);

	let code = btn.parentElement.innerHTML; // the surrounding "code" div html
	// some nasty regex to extract the code and clear tags, don't do this, kids!
	code = code.replace(/[^\S\r\n]{2,}/g, "");
	code = code.replace(/<button.*?>.*?<\/button>/gi, "");
	code = code.replace(/<span.*?>((.|\n)*?)<\/span>/gi, "$1");
	code = code.replace(/<br>(?=<br>)/gi, "\n");
	code = code.replace(/<br>/gi, "");
	code = code.replace(/&nbsp;/gi, " ");

	navigator.clipboard.writeText(code.trim());
}
