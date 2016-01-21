function scp_files($xmlInputfile,$cred){

if(!$cred){
$cred=Get-Credential
}

$doc = [XML](Get-Content $xmlInputfile);
$names=$doc.RunCollection.fileList.file.name.value;


Write-Host $name


$p=pwd;

foreach($name in $names){
Get-SCPFile -ComputerName naf-atlas.desy.de  -LocalFile "$p/$name" -RemoteFile /nfs/dust/atlas/user/peschke/analyse/$name -Credential $cred 
}

}
