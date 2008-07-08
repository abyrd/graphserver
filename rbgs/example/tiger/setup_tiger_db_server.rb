$: << "../../extension/tiger"

require 'graphserver.rb'
require 'tiger_extend.rb'

if ARGV.size < 1 then
  print "usage: ruby setup_tiger_server.rb [create_tables] [directory [directory ...] ]\n"
  exit
end

DB_PARAMS = { :host => nil,
              :port => nil,
              :options => nil,
              :tty => nil,
              :dbname => 'graphserver',
              :login => nil, #database username
              :password => nil }

gs = Graphserver.new
gs.database_params = DB_PARAMS

if ix = ARGV.index("remove_tables") then
  gs.remove_tiger_table! #clean up first

  ARGV.delete_at( ix )
end

if ix = ARGV.index("create_tables") then
  gs.create_tiger_table!

  ARGV.delete_at( ix )
end

ARGV.each do |directory|
  gs.import_tiger_to_db! directory
end