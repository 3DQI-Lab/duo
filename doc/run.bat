REM call make latexpdf
call make html

REM /y: Suppress prompt to confirm overwriting a file.
xcopy /y source\DuoManual.html build

echo %cd%
pause

