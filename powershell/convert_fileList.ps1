function convert_fileList($inputName, $thresholdFile){

$thr = Import-Csv $thresholdFile -Delimiter ",";
$doc = [XML](Get-Content $inputName);

$doc.RunCollection.fileList.file.threshold|foreach { $int_value= [convert]::ToInt64($_.value); $test= $thr[$int_value].voltage; $_.value ="$test";Write-Host $_.value }
$p=pwd;
Write-Host "$p/$inputName"
$doc.Save("$p/$inputName")

}
